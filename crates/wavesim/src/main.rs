use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand, ValueEnum};
use wavesim::{Bed, Error, RunOptions, estimate_seconds};

/// How the run is scheduled by the operating system.
#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Priority {
    /// Keep the run on the efficiency cores (macOS): quieter and cooler, but slower.
    Background,
    /// Let the system schedule it like any other program.
    Normal,
}

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
        /// Dispersive (Serre-Green-Naghdi) equations; `false` runs plain shallow water,
        /// where tall waves steepen into shocks
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        dispersive: bool,
        /// Threads for the solver, at most half the machine's cores (0 = that limit). The
        /// default of 4 keeps a laptop cool; the results are identical for any number.
        #[arg(long, default_value_t = 4)]
        threads: usize,
        /// Stop after this many minutes of wall-clock time and keep what has been simulated,
        /// as a valid shorter run (0 = no limit). A guard against a run going on too long.
        #[arg(long, default_value_t = 10.0)]
        max_minutes: f64,
        /// Scheduling priority; `background` (the default) keeps the run on macOS's
        /// efficiency cores
        #[arg(long, value_enum, default_value_t = Priority::Background)]
        priority: Priority,
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
            dispersive,
            threads,
            max_minutes,
            priority,
            out,
        } => {
            lower_priority(priority);
            let threads = use_threads(threads);
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
                dispersive,
                max_wall_seconds: (max_minutes > 0.0).then_some(max_minutes * 60.0),
            };
            if let Ok(bed) = Bed::read(&opts.bed) {
                // On the efficiency cores the same work takes about three and a half times as long
                // (measured: 15.5 s against 4.4 s for 20 simulated seconds on the Pipeline bed).
                let slower = if cfg!(target_os = "macos") && priority == Priority::Background {
                    3.5
                } else {
                    1.0
                };
                let guess = estimate_seconds(&bed, &opts, threads) * slower;
                eprintln!(
                    "estimated about {} on {threads} thread(s) (a rough guess){}",
                    human(guess),
                    match opts.max_wall_seconds {
                        Some(limit) if guess > limit => format!(
                            "; the run will probably stop at its {} limit and keep what it has",
                            human(limit)
                        ),
                        Some(limit) => format!("; it stops by itself after {}", human(limit)),
                        None => String::new(),
                    }
                );
            }
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
            if let Some(d) = summary.dispersion {
                println!(
                    "dispersive solves: {} ({:.1} iterations each, {} unconverged)",
                    d.solves,
                    d.iterations as f64 / d.solves.max(1) as f64,
                    d.unconverged
                );
            }
            if summary.truncated {
                println!("stopped at its time limit: a shorter run, but a complete one");
            }
            println!("wrote {}", summary.header.display());
            Ok(())
        }
    }
}

/// Limit the solver's worker threads, if this build has any, and say how many it will use.
///
/// Never more than half the machine's cores, so a long run cannot max out a laptop: `0`
/// means that ceiling, and a larger request is brought down to it with a note. Setting
/// `WAVESIM_ALL_CORES` lifts the ceiling for someone who really wants every core.
fn use_threads(requested: usize) -> usize {
    #[cfg(feature = "parallel")]
    {
        let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
        let ceiling = (cores / 2).max(1);
        let unlimited = std::env::var_os("WAVESIM_ALL_CORES").is_some();
        let wanted = if requested == 0 { ceiling } else { requested };
        let threads = if unlimited {
            wanted.min(cores)
        } else {
            wanted.min(ceiling)
        };
        if threads < wanted {
            eprintln!(
                "note: using {threads} threads, not {wanted}: this keeps a long run to half the \
                 machine's {cores} cores (set WAVESIM_ALL_CORES to lift the limit)"
            );
        }
        // Fails only if a pool already exists, in which case it stays as it is.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
        threads
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = requested;
        1
    }
}

/// On macOS, move this process to the background class, which the scheduler keeps on the
/// efficiency cores. Best effort: if `taskpolicy` is missing or refuses, the run goes on at
/// normal priority.
fn lower_priority(priority: Priority) {
    #[cfg(target_os = "macos")]
    if priority == Priority::Background {
        let _ = std::process::Command::new("taskpolicy")
            .args(["-b", "-p", &std::process::id().to_string()])
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status();
    }
    #[cfg(not(target_os = "macos"))]
    let _ = priority;
}

/// "45 seconds", "3.5 minutes".
fn human(seconds: f64) -> String {
    if seconds < 90.0 {
        format!("{seconds:.0} seconds")
    } else {
        format!("{:.1} minutes", seconds / 60.0)
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
