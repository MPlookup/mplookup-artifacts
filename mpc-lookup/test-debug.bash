#!/usr/bin/env bash
set -xe
trap "exit" INT TERM
trap 'kill 0' EXIT

cargo +nightly build --bin client --release

BIN=./target/release/client

n=4
max_n=4
while (( n <= max_n )); do
  echo "Running n=$n"
  "$BIN" --hosts data/2 --party 0 --mpc -n "$n" --debug 1> result_${n}_party0.txt 2>&1 & pid0=$!
  "$BIN" --hosts data/2 --party 1 --mpc -n "$n" --debug 1> result_${n}_party1.txt 2>&1 & pid1=$!
  wait "$pid0" "$pid1"
  n=$(( n * 2 ))
done

trap - INT TERM EXIT
