use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};
use wavesim::{Error, RunOptions};

#[derive(Parser)]
#[command(name = "wavesim", version, about = "Nearshore wave solver")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Run a swell over a bed file and write frames.
    Run {
        /// Bed file header, e.g. data/pipeline.json
        bed: PathBuf,
        /// Wave height, trough to crest, in metres at the wave-maker
        #[arg(long, default_value_t = 1.0)]
        height: f64,
        /// Wave period in seconds
        #[arg(long, default_value_t = 14.0)]
        period: f64,
        /// Still-water level in metres above mean sea level
        #[arg(long, default_value_t = 0.0, allow_negative_numbers = true)]
        tide: f64,
        /// Simulated seconds
        #[arg(long, default_value_t = 300.0)]
        duration: f64,
        /// Seconds between output frames
        #[arg(long, default_value_t = 2.0)]
        frame_interval: f64,
        /// Put the wave-maker where the water first gets this shallow, in metres
        #[arg(long, default_value_t = 8.0)]
        maker_depth: f64,
        /// Manning roughness in s/m^(1/3); 0 switches friction off
        #[arg(long, default_value_t = 0.0)]
        manning: f64,
        /// Output directory (default: out/<bed name>)
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Print what a bed file contains.
    Info {
        /// Bed file header
        bed: PathBuf,
    },
}

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn real_main() -> Result<(), Error> {
    match Cli::parse().command {
        Command::Info { bed } => info(&bed),
        Command::Run {
            bed,
            height,
            period,
            tide,
            duration,
            frame_interval,
            maker_depth,
            manning,
            out,
        } => {
            let name = bed
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "run".into());
            let opts = RunOptions {
                out: out.unwrap_or_else(|| PathBuf::from("out").join(name)),
                bed,
                height,
                period,
                tide,
                duration,
                frame_interval,
                maker_depth,
                manning,
            };
            let mut last_decile = -1;
            let summary = wavesim::run(&opts, |p| {
                let decile = (10.0 * p.time / p.duration) as i32;
                if decile != last_decile {
                    last_decile = decile;
                    eprintln!(
                        "  t = {:6.1} s of {:.0}  ({} frames, {} steps, highest surface {:+.2} m)",
                        p.time, p.duration, p.frames, p.steps, p.max_rise
                    );
                }
            })?;
            let plan = summary.plan;
            println!(
                "wave-maker at x = {:.0} m in {:.1} m of water (wavelength about {:.0} m); \
                 sponges {} cells offshore, {} at the sides",
                plan.maker_x,
                plan.maker_depth,
                plan.wavelength,
                plan.sponge_offshore,
                plan.sponge_side
            );
            println!(
                "{} frames, {} steps in {:.1} s; highest surface {:+.2} m above still water",
                summary.frames, summary.steps, summary.wall_seconds, summary.max_rise
            );
            println!("wrote {}", summary.header.display());
            Ok(())
        }
    }
}

fn info(path: &std::path::Path) -> Result<(), Error> {
    let bed = waveio::Bed::read(path)?;
    let h = &bed.header;
    let (lo, hi) = bed
        .elevation
        .iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(a, b), &v| {
            (a.min(v), b.max(v))
        });
    let wet =
        bed.elevation.iter().filter(|&&v| v < 0.0).count() as f64 / bed.elevation.len() as f64;
    println!("{} ({})", h.name, h.description);
    println!(
        "{} x {} cells of {} x {} m = {:.0} x {:.0} m",
        h.nx,
        h.ny,
        h.dx,
        h.dy,
        h.nx as f64 * h.dx,
        h.ny as f64 * h.dy
    );
    println!(
        "elevation {lo:.1} .. {hi:.1} m; {:.0}% below mean sea level",
        100.0 * wet
    );
    println!(
        "+x points along compass bearing {:.0} deg (the direction waves travel)",
        h.frame.x_bearing_deg
    );
    if let Some(name) = h.source.get("name").and_then(|v| v.as_str()) {
        println!("source: {name}");
    }
    Ok(())
}
