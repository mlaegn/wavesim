use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::bed::{Bed, Frame};
use crate::error::Error;

pub const RUN_FORMAT: &str = "wavesim-run";
const RUN_VERSION: u32 = 1;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RunHeader {
    pub format: String,
    pub version: u32,
    pub name: String,
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// Bed elevation file in this directory, same layout as a bed file's data.
    pub bed: String,
    /// Frames file: `frame_count` frames of every field in `fields`, back to back.
    pub frames: String,
    /// Field names in storage order within one frame. Currently `["eta"]`: free-surface
    /// elevation in metres above mean sea level (equal to the bed where it is dry).
    pub fields: Vec<String>,
    pub dtype: String,
    pub layout: String,
    pub frame_count: usize,
    /// Simulation time in seconds of each frame.
    pub times: Vec<f64>,
    pub frame: Frame,
    /// The forcing and settings the run used.
    pub waves: serde_json::Value,
    #[serde(default)]
    pub source: serde_json::Value,
}

/// Writes a run directory: `run.json`, `bed.f32` and `frames.f32`.
pub struct RunWriter {
    dir: PathBuf,
    header: RunHeader,
    frames: BufWriter<File>,
}

impl RunWriter {
    /// Create `dir`, copy the bed into it and open the frames file. `waves` records
    /// the settings so the run can be reproduced and labelled.
    pub fn create(dir: &Path, bed: &Bed, waves: serde_json::Value) -> Result<Self, Error> {
        fs::create_dir_all(dir).map_err(|e| Error::io(dir, e))?;
        let bed_path = dir.join("bed.f32");
        let bytes: Vec<u8> = bed.elevation.iter().flat_map(|v| v.to_le_bytes()).collect();
        fs::write(&bed_path, bytes).map_err(|e| Error::io(&bed_path, e))?;

        let frames_path = dir.join("frames.f32");
        let file = File::create(&frames_path).map_err(|e| Error::io(&frames_path, e))?;
        let h = &bed.header;
        Ok(Self {
            dir: dir.into(),
            header: RunHeader {
                format: RUN_FORMAT.into(),
                version: RUN_VERSION,
                name: h.name.clone(),
                nx: h.nx,
                ny: h.ny,
                dx: h.dx,
                dy: h.dy,
                bed: "bed.f32".into(),
                frames: "frames.f32".into(),
                fields: vec!["eta".into()],
                dtype: "f32-le".into(),
                layout: "per frame: row-major, y outer (increasing), x inner (increasing)".into(),
                frame_count: 0,
                times: Vec::new(),
                frame: h.frame.clone(),
                waves,
                source: h.source.clone(),
            },
            frames: BufWriter::new(file),
        })
    }

    /// Append one frame of `eta` (`nx * ny` values) taken at `time` seconds.
    pub fn write_frame(&mut self, time: f64, eta: &[f32]) -> Result<(), Error> {
        let path = self.dir.join("frames.f32");
        let want = self.header.nx * self.header.ny;
        if eta.len() != want {
            return Err(Error::format(
                &path,
                format!("frame has {} values, expected {want}", eta.len()),
            ));
        }
        for v in eta {
            self.frames
                .write_all(&v.to_le_bytes())
                .map_err(|e| Error::io(&path, e))?;
        }
        self.header.times.push(time);
        self.header.frame_count += 1;
        Ok(())
    }

    /// Flush the frames and write `run.json`. Returns the path of the header.
    pub fn finish(mut self) -> Result<PathBuf, Error> {
        let frames_path = self.dir.join("frames.f32");
        self.frames
            .flush()
            .map_err(|e| Error::io(&frames_path, e))?;
        let path = self.dir.join("run.json");
        let text = serde_json::to_string_pretty(&self.header).map_err(|source| Error::Json {
            path: path.clone(),
            source,
        })?;
        fs::write(&path, text + "\n").map_err(|e| Error::io(&path, e))?;
        Ok(path)
    }
}

/// A finished run read back from disk.
pub struct Run {
    pub header: RunHeader,
    /// All frames, concatenated: `frame_count * nx * ny` values.
    pub eta: Vec<f32>,
}

impl Run {
    pub fn read(dir: &Path) -> Result<Self, Error> {
        let path = dir.join("run.json");
        let text = fs::read_to_string(&path).map_err(|e| Error::io(&path, e))?;
        let header: RunHeader = serde_json::from_str(&text).map_err(|source| Error::Json {
            path: path.clone(),
            source,
        })?;
        if header.format != RUN_FORMAT || header.version != RUN_VERSION {
            return Err(Error::format(&path, "not a wavesim-run version 1"));
        }
        let frames_path = dir.join(&header.frames);
        let bytes = fs::read(&frames_path).map_err(|e| Error::io(&frames_path, e))?;
        let expected = header.frame_count * header.nx * header.ny * header.fields.len() * 4;
        if bytes.len() != expected {
            return Err(Error::format(
                &frames_path,
                format!("{} bytes, header promises {expected}", bytes.len()),
            ));
        }
        let eta = bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| f32::from_le_bytes(*b))
            .collect();
        Ok(Self { header, eta })
    }

    /// Frame `k` of the first field (`eta`).
    pub fn frame(&self, k: usize) -> &[f32] {
        let n = self.header.nx * self.header.ny;
        &self.eta[k * n..(k + 1) * n]
    }
}
