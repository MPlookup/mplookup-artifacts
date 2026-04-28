#!/usr/bin/env bash
set -e

# This script builds and runs a Rust binary named 'client' for benchmarking purposes.
# It iterates over increasing values of 'n' (starting from 8, doubling up to 1024),
# launching the client in multiple processes based on the number of MPC parties.
# The output from the first party (party 0) is filtered to exclude specific 'Broadcast' start/end lines,
# saved to a result file, tailed for the last 50 lines, and compressed with xz.
# Parameters:
# - parties: Controls the number of MPC nodes (defaults to 2; valid values: 1, 4, 8).
#   If parties=1, omits --mpc and --party flags, and uses data/1 for hosts.
# - test_naive: If set to 1, appends --naive flag to the command and prefixes 'naive_' to output filenames.
#   Defaults to 0 (no --naive).

parties=2  # Defaults to 2, can be set to 1, 4, or 8. Please make sure the vCPU count is 4 times parties.
test_naive=0

n=8
max_n=1024

trap "exit" INT TERM
trap 'kill 0' EXIT

cargo +nightly build --bin client --release

BIN=./target/release/client

start_pattern='^[·\s]*Start:\s+Broadcast\s+[0-9]+$'
end_pattern='^[·\s]*End:\s+Broadcast\s+[0-9]+\s+\.+[0-9]+(\.[0-9]+)?(µs|ms|ns|s)$'

mkdir -p results

while (( n <= max_n )); do
  echo "Running n=$n"
  pids=()
  for (( i=0; i<parties; i++ )); do
    cmd="$BIN --hosts data/${parties} -n $n"
    if (( parties > 1 )); then
      cmd="$cmd --party $i --mpc"
    fi
    if (( test_naive )); then
      cmd="$cmd --naive"
    fi
    if (( i == 0 )); then
      ( $cmd 2>&1 | stdbuf -o0 grep -E -v "$start_pattern|$end_pattern" | tee "results/result$( (( test_naive )) && echo "_naive" )_${parties}_parties_n_${n}.txt" ) & pids+=($!)
    else
      $cmd 1>/dev/null 2>&1 & pids+=($!)
    fi
  done
  wait "${pids[@]}"

  # tail -n 50 "results/result$( (( test_naive )) && echo "_naive" )_${parties}_parties_n_${n}.txt" > "results/result$( (( test_naive )) && echo "_naive" )_${parties}_parties_n_${n}_tail.txt"
  # xz -T0 "results/result$( (( test_naive )) && echo "_naive" )_${parties}_parties_n_${n}.txt"

  n=$(( n * 2 ))
done

trap - INT TERM EXIT