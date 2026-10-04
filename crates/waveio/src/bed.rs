use std::fs;
use std::ops::Range;
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

    /// The part of the bed in columns `x` and rows `y`, with the frame's origin moved to its
    /// lower-left corner. Panics if a range is empty or reaches past the bed.
    pub fn crop(&self, x: Range<usize>, y: Range<usize>) -> Self {
        let h = &self.header;
        assert!(
            !x.is_empty() && !y.is_empty() && x.end <= h.nx && y.end <= h.ny,
            "cannot crop columns {x:?} and rows {y:?} of a {} x {} bed",
            h.nx,
            h.ny
        );
        let elevation = y
            .clone()
            .flat_map(|j| {
                self.elevation[j * h.nx + x.start..j * h.nx + x.end]
                    .iter()
                    .copied()
            })
            .collect();
        let (along, across) = (x.start as f64 * h.dx, y.start as f64 * h.dy);
        // +x points along the bearing; +y is 90 degrees counter-clockwise from it.
        let b = h.frame.x_bearing_deg.to_radians();
        let mut header = h.clone();
        header.nx = x.len();
        header.ny = y.len();
        header.frame.origin_easting += along * b.sin() - across * b.cos();
        header.frame.origin_northing += along * b.cos() + across * b.sin();
        Self { header, elevation }
    }

    /// The bed without its first `columns` columns (the offshore end).
    pub fn crop_x(&self, columns: usize) -> Self {
        self.crop(columns..self.header.nx, 0..self.header.ny)
    }

    /// The bed on cells `factor` times as large, each the mean of the `factor` x `factor` cells
    /// it covers. Cells left over at the high-x and high-y edges are dropped; the origin stays.
    pub fn coarsen(&self, factor: usize) -> Self {
        let h = &self.header;
        let (nx, ny) = (h.nx / factor, h.ny / factor);
        assert!(
            nx > 0 && ny > 0,
            "a {} x {} bed has no {factor} x {factor} block",
            h.nx,
            h.ny
        );
        let elevation = (0..nx * ny)
            .map(|n| {
                let (i, j) = (n % nx, n / nx);
                let sum: f64 = (0..factor * factor)
                    .map(|m| {
                        let (a, b) = (i * factor + m % factor, j * factor + m / factor);
                        f64::from(self.elevation[b * h.nx + a])
                    })
                    .sum();
                (sum / (factor * factor) as f64) as f32
            })
            .collect();
        let mut header = h.clone();
        (header.nx, header.ny) = (nx, ny);
        header.dx *= factor as f64;
        header.dy *= factor as f64;
        Self { header, elevation }
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
