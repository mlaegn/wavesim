use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::Error;

pub const SLICE_FORMAT: &str = "wavesim-slice";
const SLICE_VERSION: u32 = 1;

/// A wave breaking in a vertical slice, seen side-on: the surface as a chain of nodes in every
/// frame, over a fixed bed between two walls.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SliceHeader {
    pub format: String,
    pub version: u32,
    pub name: String,
    /// The bed from the left wall to the right one, `[x, z]` in metres, `z` up from still water.
    pub bed: Vec<[f64; 2]>,
    /// Surface nodes in every frame, from the left wall to the right one.
    pub nodes: usize,
    /// Frames file: `frame_count` frames, each holding every field in `fields` for all `nodes`
    /// nodes, one field after another.
    pub frames: String,
    /// `x` and `z` are a node's position in metres; `speed` is the water's speed there in m/s.
    pub fields: Vec<String>,
    pub dtype: String,
    pub frame_count: usize,
    /// Simulation time in seconds of each frame. Frames come every few steps, and steps shorten
    /// as a lip forms, so they crowd where the wave breaks.
    pub times: Vec<f64>,
    /// The wave and settings the run used.
    pub wave: serde_json::Value,
    /// Why the run ended: `landed` (the lip reached the water), `landed against the wall...`
    /// (it did, but next to the shallow wall, which shaped it), `wall` (the wave reached the
    /// shallow wall without landing), `time limit`, or `failed: ...`.
    pub ended: String,
    /// When and where the front face first passed vertical, if it did.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub curl: Option<Curl>,
    /// When and where the lip landed, if it did.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub landing: Option<Landing>,
}

/// The moment the front face first passes vertical.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Curl {
    pub time: f64,
    /// Where the face is vertical, m.
    pub x: f64,
    /// The crest's height above still water then, m.
    pub crest: f64,
    /// The still-water depth under the vertical face, m.
    pub depth: f64,
}

/// The moment the lip reaches the water below it.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Landing {
    pub time: f64,
    /// Where the lip's tip lands, m.
    pub x: f64,
    pub z: f64,
    /// How far forward of where the face first stood vertical the lip lands, m.
    pub throw: f64,
    /// The air the lip closes off, m^2 per metre of crest.
    pub tube_area: f64,
    /// The tube's reach from back to front and from bottom to top, m.
    pub tube_width: f64,
    pub tube_height: f64,
}

/// Writes a slice directory: `slice.json` and `frames.f32`.
pub struct SliceWriter {
    dir: PathBuf,
    header: SliceHeader,
    frames: BufWriter<File>,
}

impl SliceWriter {
    /// Create `dir` and open the frames file, for `nodes` surface nodes over `bed`.
    pub fn create(
        dir: &Path,
        bed: &[(f64, f64)],
        nodes: usize,
        wave: serde_json::Value,
    ) -> Result<Self, Error> {
        fs::create_dir_all(dir).map_err(|e| Error::io(dir, e))?;
        let path = dir.join("frames.f32");
        let file = File::create(&path).map_err(|e| Error::io(&path, e))?;
        let name = dir
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "slice".into());
        let header = SliceHeader {
            format: SLICE_FORMAT.into(),
            version: SLICE_VERSION,
            name,
            bed: bed.iter().map(|p| [p.0, p.1]).collect(),
            nodes,
            frames: "frames.f32".into(),
            fields: vec!["x".into(), "z".into(), "speed".into()],
            dtype: "<f4".into(),
            frame_count: 0,
            times: Vec::new(),
            wave,
            ended: "unfinished".into(),
            curl: None,
            landing: None,
        };
        Ok(Self {
            dir: dir.to_path_buf(),
            header,
            frames: BufWriter::new(file),
        })
    }

    /// Append a frame: each node's position and the water's speed there.
    pub fn write_frame(
        &mut self,
        time: f64,
        surface: &[(f64, f64)],
        speed: &[f64],
    ) -> Result<(), Error> {
        let path = self.dir.join("frames.f32");
        if surface.len() != self.header.nodes || speed.len() != self.header.nodes {
            return Err(Error::format(
                &path,
                format!(
                    "a frame of {} nodes, expected {}",
                    surface.len(),
                    self.header.nodes
                ),
            ));
        }
        let values = surface
            .iter()
            .map(|p| p.0)
            .chain(surface.iter().map(|p| p.1))
            .chain(speed.iter().copied());
        for v in values {
            self.frames
                .write_all(&(v as f32).to_le_bytes())
                .map_err(|e| Error::io(&path, e))?;
        }
        self.header.frame_count += 1;
        self.header.times.push(time);
        Ok(())
    }

    /// Frames written so far.
    pub fn frame_count(&self) -> usize {
        self.header.frame_count
    }

    /// Flush the frames and write `slice.json`, saying how the run ended and what happened.
    pub fn finish(
        mut self,
        ended: &str,
        curl: Option<Curl>,
        landing: Option<Landing>,
    ) -> Result<PathBuf, Error> {
        let frames_path = self.dir.join("frames.f32");
        self.frames
            .flush()
            .map_err(|e| Error::io(&frames_path, e))?;
        self.header.ended = ended.into();
        self.header.curl = curl;
        self.header.landing = landing;
        let path = self.dir.join("slice.json");
        let text = serde_json::to_string_pretty(&self.header).map_err(|source| Error::Json {
            path: path.clone(),
            source,
        })?;
        fs::write(&path, text + "\n").map_err(|e| Error::io(&path, e))?;
        Ok(path)
    }
}

/// A slice read back: its header and every frame's values.
pub struct Slice {
    pub header: SliceHeader,
    pub data: Vec<f32>,
}

impl Slice {
    pub fn read(dir: &Path) -> Result<Self, Error> {
        let path = dir.join("slice.json");
        let text = fs::read_to_string(&path).map_err(|e| Error::io(&path, e))?;
        let header: SliceHeader = serde_json::from_str(&text).map_err(|source| Error::Json {
            path: path.clone(),
            source,
        })?;
        if header.format != SLICE_FORMAT || header.version != SLICE_VERSION {
            return Err(Error::format(&path, "not a wavesim-slice version 1"));
        }
        if header.times.len() != header.frame_count {
            return Err(Error::format(&path, "one time per frame"));
        }
        let frames_path = dir.join(&header.frames);
        let bytes = fs::read(&frames_path).map_err(|e| Error::io(&frames_path, e))?;
        let expected = header.frame_count * header.nodes * header.fields.len() * 4;
        if bytes.len() != expected {
            return Err(Error::format(
                &frames_path,
                format!("{} bytes, header promises {expected}", bytes.len()),
            ));
        }
        let data = bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| f32::from_le_bytes(*b))
            .collect();
        Ok(Self { header, data })
    }

    /// Field `field` (an index into `fields`) of frame `frame`, one value per node.
    pub fn field(&self, frame: usize, field: usize) -> &[f32] {
        let n = self.header.nodes;
        let start = (frame * self.header.fields.len() + field) * n;
        &self.data[start..start + n]
    }
}
