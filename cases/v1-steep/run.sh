#!/bin/bash
# Run V1 (docs/3d-box.md) at one cell size, guarded so the laptop stays quiet:
#
#   cases/v1-steep/run.sh 0.4 1200     # cell size in metres, wall-clock limit in seconds
#   cases/v1-steep/run.sh 0.1 5400 --iso --from 0.1@14.5        # interIsoFoam from that run
#   cases/v1-steep/run.sh 0.05 10800 --x0 170 --from 0.1@14.5   # finer, the tank cut at 170 m
#
# OpenFOAM runs in Docker on one core. Every 30 s a watchdog reads macOS's thermal state and
# stops the run at anything above nominal, or once the limit has passed; caffeinate keeps the
# Mac awake only while this script runs. The case goes to out/foam-v1-steep-<dx>/, its log of
# thermal readings to watchdog.log there; a variant's case is named after its flags too
# (out/foam-v1-steep-0.1-iso/, out/foam-v1-steep-0.05-x170/). Exit status: 0 finished, 2 stopped by heat, 3 by time.
set -euo pipefail

dx=${1:?cell size in metres}
limit=${2:?wall-clock limit in seconds}
shift 2
solver=interFoam suffix="" mesh_args=() from=""
while [ $# -gt 0 ]; do
    case $1 in
        --iso) solver=interIsoFoam; suffix+="-iso"; shift ;;
        --x0) mesh_args+=(--x0 "$2"); suffix+="-x$2"; shift 2 ;;
        --from) from=$2; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
done
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
case_dir="$root/out/foam-v1-steep-$dx$suffix"
image=opencfd/openfoam-default:2512
name="v1-steep-${dx//./_}$suffix"
thermal="$root/out/foam-tools/thermal"

caffeinate -is -w $$ &

if [ ! -x "$thermal" ]; then
    mkdir -p "$(dirname "$thermal")"
    swiftc -O "$here/thermal.swift" -o "$thermal"
fi
docker info >/dev/null 2>&1 || {
    open -a OrbStack
    for _ in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 2; done
}
if docker ps -q --filter "name=^${name}$" | grep -q .; then
    echo "a container named $name is already running" >&2
    exit 1
fi

py() { (cd "$here" && uv run -q --python 3.12 --with numpy python -I make_case.py "$@"); }
# The image's login shell starts in its home folder, so every command cds to the case first.
foam() { docker run --rm --cpus 1 -m 4g -v "$case_dir:/case" "$image" bash -lc "cd /case && $1"; }

if [ ! -f "$case_dir/0/C" ]; then
    rm -rf "$case_dir"
    py mesh --dx "$dx" --solver "$solver" ${mesh_args[@]+"${mesh_args[@]}"} --out "$case_dir"
    foam 'blockMesh > log.blockMesh 2>&1 && checkMesh > log.checkMesh 2>&1 &&
          postProcess -func writeCellCentres -time 0 > log.centres 2>&1'
    if [ -n "$from" ]; then
        py restart --from "$root/out/foam-v1-steep-${from%@*}" --time "${from#*@}" --out "$case_dir" \
            | tee "$case_dir/log.init"
    else
        py init --out "$case_dir" | tee "$case_dir/log.init"
    fi
fi

log="$case_dir/watchdog.log"
docker run -d --rm --name "$name" --cpus 1 -m 4g -v "$case_dir:/case" "$image" \
    bash -lc "cd /case && $solver > log.interFoam 2>&1" >/dev/null
start=$(date +%s)
echo "$(date '+%F %T') started $name, limit ${limit}s" | tee -a "$log"
status=0
while docker ps -q --filter "name=^${name}$" | grep -q .; do
    state=$("$thermal" 2>/dev/null || echo unreadable)
    elapsed=$(( $(date +%s) - start ))
    # || true: the log does not exist for the first moments, and set -e must not end the watch.
    sim=$( (grep -a '^Time = ' "$case_dir/log.interFoam" 2>/dev/null || true) | tail -1 | cut -d' ' -f3)
    echo "$(date '+%F %T') ${elapsed}s thermal=$state sim_time=${sim:-0}" >> "$log"
    if [ "$state" != nominal ]; then
        echo "$(date '+%F %T') thermal state $state: stopping" | tee -a "$log"
        docker stop -t 5 "$name" >/dev/null; status=2; break
    fi
    if [ "$elapsed" -gt "$limit" ]; then
        echo "$(date '+%F %T') over ${limit}s: stopping" | tee -a "$log"
        docker stop -t 5 "$name" >/dev/null; status=3; break
    fi
    sleep 30
done
echo "$(date '+%F %T') ended after $(( $(date +%s) - start ))s, status $status; containers left: $(docker ps -q | wc -l | tr -d ' ')" | tee -a "$log"
exit $status
