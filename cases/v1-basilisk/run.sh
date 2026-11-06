#!/bin/bash
# Run V1 (docs/3d-box.md) in Basilisk at one refinement level, guarded as cases/v1-steep/run.sh
# is, so the laptop stays quiet:
#
#   cases/v1-basilisk/run.sh 11 3600        # level (258 m / 2^level cells), wall-clock limit in s
#   cases/v1-basilisk/run.sh 10 3600 --3d   # the 3D strip of v1-3d.c, from the 2D run's state
#
# Basilisk must be built in out/basilisk/basilisk/src (see v1.c). The case is compiled and run
# in Docker on one core, in out/basilisk-v1-L<level>/ (out/basilisk-v1-3d-L<level>/ for --3d,
# which needs the 2D run at that level first). Every 30 s a watchdog reads macOS's
# thermal state and stops the run at anything above nominal, or once the limit has passed;
# caffeinate keeps the Mac awake only while this script runs. Exit status: 0 finished,
# 2 stopped by heat, 3 by time.
set -euo pipefail

level=${1:?refinement level}
limit=${2:?wall-clock limit in seconds}
dims=${3:-}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
src=v1 tag=""
[ "$dims" = --3d ] && src=v1-3d tag=-3d
case_dir="$root/out/basilisk-v1$tag-L$level"
image=opencfd/openfoam-default:2512
name="basilisk-v1$tag-L$level"
thermal="$root/out/foam-tools/thermal"

caffeinate -is -w $$ &

if [ ! -x "$thermal" ]; then
    mkdir -p "$(dirname "$thermal")"
    swiftc -O "$root/cases/v1-steep/thermal.swift" -o "$thermal"
fi
docker info >/dev/null 2>&1 || {
    open -a OrbStack
    for _ in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 2; done
}
if docker ps -q --filter "name=^${name}$" | grep -q .; then
    echo "a container named $name is already running" >&2
    exit 1
fi

rm -rf "$case_dir"
mkdir -p "$case_dir"
cp "$here/$src.c" "$case_dir/"
if [ -n "$tag" ]; then
    cp "$root/out/basilisk-v1-L$level/state-14.500.bin" "$case_dir/state.bin"
fi
# Basilisk's qcc finds its headers where it was built, so the build is mounted at that path.
docker run -d --rm --name "$name" --cpus 1 -m 4g \
    -v "$root/out/basilisk:/bas:ro" -v "$case_dir:/case" "$image" \
    bash -c "export BASILISK=/bas/basilisk/src PATH=/bas/basilisk/src:\$PATH; cd /case &&
             qcc -O2 -Wall -disable-dimensions $src.c -o $src -lm > log.build 2>&1 &&
             ./$src $level > log.out 2> log.err" >/dev/null
log="$case_dir/watchdog.log"
start=$(date +%s)
echo "$(date '+%F %T') started $name, limit ${limit}s" | tee -a "$log"
status=0
while docker ps -q --filter "name=^${name}$" | grep -q .; do
    state=$("$thermal" 2>/dev/null || echo unreadable)
    elapsed=$(( $(date +%s) - start ))
    sim=$( (tail -1 "$case_dir/log" 2>/dev/null || true) | cut -d' ' -f2)
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
