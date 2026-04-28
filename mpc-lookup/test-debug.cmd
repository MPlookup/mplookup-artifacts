del target\debug\client.exe
cargo +nightly build --bin client
set RUST_BACKTRACE=full
start "" cmd /k "target\debug\client --hosts data/2 --party 0 --mpc -n 4 --debug"
start "" cmd /k "target\debug\client --hosts data/2 --party 1 --mpc -n 4 --debug"

pause