#!/bin/sh
# Every test, including the slow physics checks that run whole simulations, without heating the
# laptop: on macOS on the efficiency cores (background priority), one test at a time, with two
# solver threads, compiling on two jobs. It takes a while; start it and leave it.
#
#     scripts/slow-tests.sh            # everything
#     scripts/slow-tests.sh lab        # only tests whose names match "lab"
set -eu
cd "$(dirname "$0")/.."
export RUST_TEST_THREADS=1 RAYON_NUM_THREADS=2 CARGO_BUILD_JOBS=2
if command -v taskpolicy >/dev/null 2>&1; then
    exec taskpolicy -b cargo test --release -- --include-ignored "$@"
else
    exec nice -n 19 cargo test --release -- --include-ignored "$@"
fi
