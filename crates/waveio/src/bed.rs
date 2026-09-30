use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};
use wavecore::{Bathymetry, Grid};

use crate::error::Error;

pub const BED_FORMAT: &str = "wavesim-bed";
const BED_VERSION: u32 = 1;
const DTYPE: &str = "f32-le";

/// Where the grid sits on the Earth and how it is rotated.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    pub crs: String,
    pub origin_easting: f64,
    pub origin_northing: f64,
    /// Compass bearing (degrees clockwise from north) that +x points along.
    pub x_bearing_deg: f64,
    #[serde(default)]
    pub note: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BedHeader {
    pub format: String,
    pub version: u32,
    pub name: String,
    #[serde(default)]
    pub description: String,
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// File name of the raw data, relative to the header.
    pub data: String,
    pub dtype: String,
    pub layout: String,
    pub units: String,
    pub frame: Frame,
    #[serde(default)]
    pub stats: serde_json::Value,
    #[serde(default)]
    pub source: serde_json::Value,
}

/// A bed on disk: header plus bed elevation in metres above mean sea level.
#[derive(Clone, Debug)]
pub struct Bed {
    pub header: BedHeader,
    /// Row-major, `y` outer and increasing, `x` inner and increasing.
    pub elevation: Vec<f32>,
}

impl Bed {
    /// Read `header_path` and the data file it names. Rejects anything that is not a
    /// well-formed version-1 bed: wrong format, wrong length, NaN or infinite values.
    pub fn read(header_path: &Path) -> Result<Self, Error> {
        let text = fs::read_to_string(header_path).map_err(|e| Error::io(header_path, e))?;
        let header: BedHeader = serde_json::from_str(&text).map_err(|source| Error::Json {
            path: header_path.into(),
            source,
        })?;
        if header.format != BED_FORMAT || header.version != BED_VERSION {
            return Err(Error::format(
                header_path,
                format!(
                    "expected {BED_FORMAT} version {BED_VERSION}, found {} version {}",
                    header.format, header.version
                ),
            ));
        }
        if header.dtype != DTYPE {
            return Err(Error::format(
                header_path,
                format!("unsupported dtype {}, expected {DTYPE}", header.dtype),
            ));
        }
        if header.nx == 0 || header.ny == 0 || header.dx <= 0.0 || header.dy <= 0.0 {
            return Err(Error::format(
                header_path,
                "grid size and spacing must be positive",
            ));
        }

        let data_path = header_path.with_file_name(&header.data);
        let bytes = fs::read(&data_path).map_err(|e| Error::io(&data_path, e))?;
        let expected = header.nx * header.ny * 4;
        if bytes.len() != expected {
            return Err(Error::format(
                &data_path,
                format!(
                    "{} bytes, but {} x {} f32 values need {expected}",
                    bytes.len(),
                    header.nx,
                    header.ny
                ),
            ));
        }
        let elevation: Vec<f32> = bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| f32::from_le_bytes(*b))
            .collect();
        if let Some(n) = elevation.iter().position(|v| !v.is_finite()) {
            return Err(Error::format(
                &data_path,
                format!("value {n} is not finite; the bed must have no gaps"),
            ));
        }
        Ok(Self { header, elevation })
    }

    /// Write the header to `header_path` and the data next to it under `header.data`.
    pub fn write(&self, header_path: &Path) -> Result<(), Error> {
        let data_path = header_path.with_file_name(&self.header.data);
        if let Some(dir) = header_path.parent() {
            fs::create_dir_all(dir).map_err(|e| Error::io(dir, e))?;
        }
        let bytes: Vec<u8> = self
            .elevation
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        fs::write(&data_path, bytes).map_err(|e| Error::io(&data_path, e))?;
        let text = serde_json::to_string_pretty(&self.header).map_err(|source| Error::Json {
            path: header_path.into(),
            source,
        })?;
        fs::write(header_path, text + "\n").map_err(|e| Error::io(header_path, e))
    }

    /// A version-1 header for a bed with the given geometry; the frame is a placeholder
    /// for beds that do not come from real data.
    pub fn header_for(name: &str, nx: usize, ny: usize, dx: f64, dy: f64) -> BedHeader {
        BedHeader {
            format: BED_FORMAT.into(),
            version: BED_VERSION,
            name: name.into(),
            description: String::new(),
            nx,
            ny,
            dx,
            dy,
            data: format!("{name}.f32"),
            dtype: DTYPE.into(),
            layout: "row-major, y outer (increasing), x inner (increasing)".into(),
            units: "metres above mean sea level; negative is underwater".into(),
            frame: Frame {
                crs: "local".into(),
                origin_easting: 0.0,
                origin_northing: 0.0,
                x_bearing_deg: 0.0,
                note: String::new(),
            },
            stats: serde_json::Value::Null,
            source: serde_json::Value::Null,
        }
    }

    pub fn grid(&self) -> Grid {
        Grid::new(
            self.header.nx,
            self.header.ny,
            self.header.dx,
            self.header.dy,
        )
    }

    /// The bed as the solver sees it, with the still-water level moved to zero: pass the
    /// tide as `still_water` (metres above mean sea level) to run at that tide.
    pub fn bathymetry(&self, still_water: f64) -> Bathymetry {
        let values: Vec<f64> = self
            .elevation
            .iter()
            .map(|&v| f64::from(v) - still_water)
            .collect();
        Bathymetry::from_interior(&self.grid(), &values)
    }
}
