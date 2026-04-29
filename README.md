# Artifacts for 'MPlookup'

This repository contains artifacts for the paper 'MPlookup'. These artifacts are currently provided anonymously, exclusively for peer review purposes. The authors, who remain anonymous due to submission guidelines, will publicly release these artifacts with their identities following the completion of the review process.

We developed the MPlookup codebase based on [collaborative-zksnark](https://github.com/alex-ozdemir/collaborative-zksnark). We also integrated a Rust implementation of MPC protocols from the C# [CompatCircuit](https://github.com/BDS-SDU/vdoram-artifacts) implementation.

## LICENSE

The license file for the anonymous artifacts.

## mpc-lookup

The main component of the MPlookup system. It includes:

- Our implementation of $O(N \log^2 N)$ MPlookup. It also includes the $O(N^2)$ strawman for comparison. Location: `mpc-lookup/src`.
- A simplified, easy-to-read single-party emulation of the MPlookup preprocessing algorithm. Location: `mpc-lookup/references/oblivious_lookup_permutation.single-party-emulation.py`.
- A script to run the experiments. Location: `mpc-lookup/test.bash`.
- Collected experiment results. Location: `mpc-lookup/results`.
- A script to draw the experiment figures in the Evaluation section. Location: `mpc-lookup/result-analyze/draw.py`.
- A script to produce all the numerical data used in the Evaluation section. Location: `mpc-lookup/result-analyze/compute-evaluation-data.py`.

Instructions to run the experiments:

### 1. System Requirements

A machine with the following minimum specifications is required:

- Operating System: Debian 13 or Ubuntu Server 24.04 LTS (amd64)
- CPU: 8 cores (4 cores per party)
- Memory: 16 GB

### 2. Rust Nightly Toolchain Installation

The project requires the nightly toolchain of the Rust programming language. The following command installs rustup (the Rust toolchain manager) and configures the nightly version as the default only if there are no existing Rust installations.

```bash
curl --proto '=https' --tlsv1.3 https://sh.rustup.rs -sSf | sh -s -- --default-toolchain nightly -y
```

Note for users in Mainland China: To utilize a regional mirror for the Rust installation, execute the following commands prior to the curl command above:

```bash
echo 'export RUSTUP_UPDATE_ROOT=https://mirrors.tuna.tsinghua.edu.cn/rustup/rustup' >> ~/.bashrc
echo 'export RUSTUP_DIST_SERVER=https://mirrors.tuna.tsinghua.edu.cn/rustup' >> ~/.bashrc
source ~/.bashrc

mkdir -vp ${CARGO_HOME:-$HOME/.cargo}

cat << EOF | tee -a ${CARGO_HOME:-$HOME/.cargo}/config.toml
[source.crates-io]
replace-with = 'mirror'

[source.mirror]
registry = "sparse+https://mirrors.tuna.tsinghua.edu.cn/crates.io-index/"

[registries.mirror]
index = "sparse+https://mirrors.tuna.tsinghua.edu.cn/crates.io-index/"
EOF
```

A successful installation will produce the following output:

```
info: default toolchain set to 'nightly-x86_64-unknown-linux-gnu'

  nightly-x86_64-unknown-linux-gnu installed - rustc 1.90.0-nightly (da58c0513 2025-07-03)


Rust is installed now. Great!
```

### 3. Run the experiments

Edit the experiment script `test.bash`:

```sh
cd mpc-lookup/
nano test.bash
```

The content of `test.bash`:

```sh
parties=2  # Defaults to 2, can be set to 1, 2, 4, 8, 16. Please make sure the vCPU count is 4 times the number of parties.
test_naive=0

n=8
max_n=1024
```

Below explains the parameters:
- `parties`: the number of multi-party computation parties, i.e., provers.
- `test_naive`: whether to test the naive implementation (0: test MPlookup; 1: test the strawman).
- `n`: the initial size of the input vector t. The script will double the size of t in each iteration until it reaches `max_n`. Must be a power of 2.
- `max_n`: the maximum size of the input vector t. Must be a power of 2.

Then run the script. We suggest using `tmux` to run the script, as it may take a long time to finish.

```sh
bash test.bash
```

In the end, the script will save the log files to the `mpc-lookup/results` folder. Each log file contains the running time of each step of the MPlookup protocol, which can be used to draw the experiment figures in the Evaluation section.

### 4. Run the evaluation data generation and figure drawing scripts

To generate the evaluation data or draw figures used in the Evaluation section, run the following commands:

Initially, set up a Python virtual environment and install the required dependencies:

```sh
cd mpc-lookup/result-analyze/
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Edit the figure format in `draw.py` if needed. By default, it will save the figures in PDF format.

Then, run the scripts to compute evaluation data and draw figures:

```sh
python compute-evaluation-data.py
python draw.py
```
