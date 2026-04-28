// Example usage:
// cargo +nightly run --bin client -- --hosts hosts.txt --party 0 --mpc -n 128
// cargo +nightly run --bin client -- -n 128

use ark_std::{end_timer, start_timer, test_rng, UniformRand};
use rand::Rng;
use std::path::PathBuf;
use structopt::StructOpt;

use ark_bls12_377::Fr as F;
use mpc_algebra::honest_but_curious::{MpcField, MpcPairingEngine};
use mpc_algebra::Reveal;
use mpc_net::{MpcMultiNet, MpcNet};

use mpc_lookup::{secure_oblivious_lookup_permutation, naive_secure_oblivious_lookup_permutation, verify_lookup_permutation_outputs_in_plaintext};

use ark_bls12_377::Bls12_377;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly_commit::marlin::marlin_pc::{CommitterKey, MarlinKZG10};
use ark_poly_commit::PolynomialCommitment;
use blake2::Blake2s;
use mpc_lookup::halo2_lookup::{
    reveal_proof, LookupConfig, LookupOpenings, LookupProof, LookupProver,
    LookupTranscript, LookupVerificationKey, LookupVerifier,
};

type MF = MpcField<F>;
type ME = MpcPairingEngine<Bls12_377>;
type LocalPC = MarlinKZG10<Bls12_377, DensePolynomial<F>>;

#[derive(StructOpt, Debug)]
struct Opt {
    /// Hosts list file (for MPC mode)
    #[structopt(long, parse(from_os_str))]
    hosts: Option<PathBuf>,

    /// Our party index (for MPC mode)
    #[structopt(long, default_value = "0")]
    party: usize,

    /// Use MPC (otherwise run the public/local variant)
    #[structopt(long)]
    mpc: bool,

    /// Table size n
    #[structopt(short, default_value = "4")]
    n: usize,

    /// Number of queries m (default n/2)
    #[structopt(short)]
    m: Option<usize>,

    /// Enable debug output (reveals MPC values)
    #[structopt(long)]
    debug: bool,

    /// Use naive permutation method instead of efficient method
    #[structopt(long)]
    naive: bool,
}

fn prepare_input<R: Rng>(n: usize, m: usize, rng: &mut R) -> (Vec<F>, Vec<F>) {
    let mut t: Vec<F> = Vec::with_capacity(n);
    while t.len() < n {
        let x = F::rand(rng);
        if !t.iter().any(|&y| y == x) {
            t.push(x);
        }
    }
    let mut f: Vec<F> = (0..m).map(|_| t[(rng.next_u64() as usize) % n]).collect();
    {
        let mut tmp = f.clone();
        tmp.sort();
        let has_dup = tmp.windows(2).any(|w| w[0] == w[1]);
        if !has_dup && m >= 2 {
            f[m - 1] = f[0];
        }
    }
    (t, f)
}

fn pad_f_prime_to_match_t_prime(
    mut f_prime: Vec<MF>,
    t_prime: &[MF],
    debug: bool,
) -> Vec<MF> {
    if f_prime.len() < t_prime.len() {
        if debug {
            println!(
                "Padding f' from length {} to {} using corresponding t' values",
                f_prime.len(),
                t_prime.len()
            );
        }
        for i in f_prime.len()..t_prime.len() {
            f_prime.push(t_prime[i].clone());
        }
    }
    f_prime
}

fn pad_f_to_match_t(
    mut f: Vec<MF>,
    t: &[MF],
    t_prime: &[MF],
    debug: bool,
) -> Vec<MF> {
    if f.len() < t.len() {
        if debug {
            println!(
                "Padding A (f) from length {} to {} using corresponding t' values",
                f.len(),
                t.len()
            );
        }
        for i in f.len()..t.len() {
            f.push(t_prime[i].clone());
        }
    }
    f
}

fn setup_halo2_proof_system(
    n: usize,
    m: usize,
    debug: bool,
) -> (LookupConfig, CommitterKey<ME>, ark_poly_commit::marlin::marlin_pc::VerifierKey<Bls12_377>, std::time::Duration) {
    let setup_start = std::time::Instant::now();
    let setup_timer = start_timer!(|| "MarlinKZG10 setup");
    let config = LookupConfig::new(n, m, 5);
    let mut setup_rng = test_rng();
    let max_degree = config.domain_size * 3;
    let pp = LocalPC::setup(max_degree, None, &mut setup_rng)
        .expect("Failed to setup MarlinKZG10");
    let (ck_local, vk_pcs_local) =
        LocalPC::trim(&pp, max_degree, 0, Some(&[config.domain_size - 1]))
            .expect("Failed to trim MarlinKZG10 parameters");
    if debug {
        eprintln!("=== SETUP DEBUG ===");
        eprintln!("max_degree (universal params): {}", max_degree);
        eprintln!("domain_size: {}", config.domain_size);
        eprintln!("supported_degrees: [{}]", config.domain_size - 1);
    }
    let ck = <CommitterKey<ME> as Reveal>::from_public(ck_local);
    end_timer!(setup_timer);
    let setup_time = setup_start.elapsed();
    (config, ck, vk_pcs_local, setup_time)
}

fn generate_halo2_proof(
    config: LookupConfig,
    ck: CommitterKey<ME>,
    a_values: Vec<MF>,
    s_values: Vec<MF>,
    t_prime_values: Vec<MF>,
    f_prime_values: Vec<MF>,
    debug: bool,
) -> (
    LookupProof<
        F,
        <LocalPC as PolynomialCommitment<F, DensePolynomial<F>>>::Commitment,
        <LocalPC as PolynomialCommitment<F, DensePolynomial<F>>>::Proof,
    >,
    std::time::Duration,
) {
    let prove_start = std::time::Instant::now();
    let prove_timer = start_timer!(|| "Proof generation and reveal");
    let mut transcript_p = LookupTranscript::<Blake2s>::new(b"halo2-lookup");
    let mut prover = LookupProver::<Blake2s>::new(config, ck, a_values, s_values, debug);
    
    prover.set_permuted_values(t_prime_values, f_prime_values);
    let proof = prover.prove(&mut transcript_p);
    
    println!("Proof generated successfully");
    
    let revealed_proof = LookupProof {
        commitments: proof.commitments,
        openings: LookupOpenings {
            a_zeta: proof.openings.a_zeta,
            s_zeta: proof.openings.s_zeta,
            a_prime_zeta: proof.openings.a_prime_zeta,
            s_prime_zeta: proof.openings.s_prime_zeta,
            z_zeta: proof.openings.z_zeta,
            z_omega_zeta: proof.openings.z_omega_zeta,
            a_prime_omega_inv_zeta: proof.openings.a_prime_omega_inv_zeta,
            q_blind_zeta: proof.openings.q_blind_zeta,
            q_last_zeta: proof.openings.q_last_zeta,
            l_0_zeta: proof.openings.l_0_zeta,
            q_zeta: proof.openings.q_zeta,
            proof_zeta: reveal_proof(proof.openings.proof_zeta),
            proof_omega_zeta: reveal_proof(proof.openings.proof_omega_zeta),
            proof_omega_inv_zeta: reveal_proof(proof.openings.proof_omega_inv_zeta),
        },
    };
    end_timer!(prove_timer);
    let prove_time = prove_start.elapsed();
    (revealed_proof, prove_time)
}

fn verify_halo2_proof(
    config: LookupConfig,
    vk_pcs_local: ark_poly_commit::marlin::marlin_pc::VerifierKey<Bls12_377>,
    proof: &LookupProof<
        F,
        <LocalPC as PolynomialCommitment<F, DensePolynomial<F>>>::Commitment,
        <LocalPC as PolynomialCommitment<F, DensePolynomial<F>>>::Proof,
    >,
    debug: bool,
) -> (bool, std::time::Duration) {
    let verify_start = std::time::Instant::now();
    let verify_timer = start_timer!(|| "Proof verification");
    let mut transcript_v = LookupTranscript::<Blake2s>::new(b"halo2-lookup");
    let domain = ark_poly::domain::radix2::Radix2EvaluationDomain::<F>::new(config.domain_size)
        .expect("Failed to create domain");
    let vk = LookupVerificationKey {
        omega: domain.group_gen,
        domain_size: config.domain_size,
        com_q_blind: None,
        com_q_last: None,
        config,
    };
    let verifier = LookupVerifier::<F, DensePolynomial<F>, LocalPC, Blake2s>::new(vk, vk_pcs_local, debug);
    let is_valid = verifier.verify(proof, &mut transcript_v);
    end_timer!(verify_timer);
    let verify_time = verify_start.elapsed();
    
    if debug {
        println!(
            "Proof verification result: {}",
            if is_valid { "VALID" } else { "INVALID" }
        );
    }
    
    (is_valid, verify_time)
}

fn main() {
    let opt = Opt::from_args();
    let n = opt.n;
    let m = opt.m.unwrap_or(n / 2);

    if opt.mpc {
        let hosts_path = opt
            .hosts
            .as_ref()
            .expect("MPC mode requires --hosts path to host file");
        MpcMultiNet::init_from_file(hosts_path.to_str().unwrap(), opt.party);
        println!("Initialized MPC network from {:?}", hosts_path);
    } else {
        println!("Running local/public variant");
    }

    if opt.mpc {
        let mut rng = test_rng();
        let (t, f) = prepare_input(n, m, &mut rng);
        let t_plain = t.clone();
        let f_plain = f.clone();

        let t: Vec<MF> = MF::king_share_batch(t, &mut rng);
        let f: Vec<MF> = MF::king_share_batch(f, &mut rng);

        let method_name = if opt.naive {
            "naive_secure_oblivious_lookup_permutation"
        } else {
            "secure_oblivious_lookup_permutation"
        };

        println!(
            "Running {} in MPC mode (n={}, m={})",
            method_name, n, m
        );

        let permutation_start_time = std::time::Instant::now();
        let permutation_timer =
            start_timer!(|| format!("{} MPC n={}, m={}", method_name, n, m));
        let (t_prime_mpc, f_prime_mpc) = if opt.naive {
            naive_secure_oblivious_lookup_permutation(n, m, &t, &f, opt.debug)
        } else {
            secure_oblivious_lookup_permutation(n, m, &t, &f, opt.debug)
        };
        end_timer!(permutation_timer);
        let permutation_time_elapsed = permutation_start_time.elapsed();

        let f_prime_mpc = pad_f_prime_to_match_t_prime(f_prime_mpc, &t_prime_mpc, opt.debug);

        if opt.debug {
            let t_prime_pub: Vec<F> = t_prime_mpc.iter().map(|x| x.reveal()).collect();
            let f_prime_pub_unpadded: Vec<F> = f_prime_mpc.iter().take(m).map(|x| x.reveal()).collect();
            let f_prime_pub: Vec<F> = f_prime_mpc.iter().map(|x| x.reveal()).collect();

            println!(
                "t' length = {}, f' length = {} (padded)",
                t_prime_pub.len(),
                f_prime_pub.len()
            );

            if n <= 16 {
                println!("\nDEBUG: After {}:", method_name);
                println!("f (original queries) = {:?}", f_plain.iter().take(8).collect::<Vec<_>>());
                println!("f' (sorted queries, unpadded) = {:?}", f_prime_pub_unpadded.iter().take(8).collect::<Vec<_>>());
                println!("f' (sorted queries, padded) = {:?}", f_prime_pub.iter().take(8).collect::<Vec<_>>());
                println!("t (original table) = {:?}", t_plain.iter().take(8).collect::<Vec<_>>());
                println!("t' (permuted table) = {:?}", t_prime_pub.iter().take(8).collect::<Vec<_>>());
            }

            verify_lookup_permutation_outputs_in_plaintext(n, m, &t_plain, &f_plain, &t_prime_pub, &f_prime_pub_unpadded, opt.debug);

            if n <= 16 {
                for (i, v) in t_prime_pub.iter().enumerate() {
                    println!("t'[{}] = {}", i, v);
                }
                for (i, v) in f_prime_pub.iter().enumerate() {
                    println!("f'[{}] = {}", i, v);
                }
            }
        }

        let (setup_time, prove_time, verify_time) = {
            println!("\n=== Halo2 Lookup Argument PROOF Generation (MPC-Compatible) ===");
            
            let (config, ck, vk_pcs_local, setup_time) = setup_halo2_proof_system(n, m, opt.debug);
            
            let f_padded = pad_f_to_match_t(f.clone(), &t, &t_prime_mpc, opt.debug);
            let s_values = t.clone();
            let a_values = f_padded;
            
            let (revealed_proof, prove_time) = generate_halo2_proof(
                config.clone(),
                ck,
                a_values,
                s_values,
                t_prime_mpc.clone(),
                f_prime_mpc.clone(),
                opt.debug,
            );
            
            let (is_valid, verify_time) = verify_halo2_proof(config, vk_pcs_local, &revealed_proof, opt.debug);
            
            if !is_valid {
                eprintln!("WARNING: Proof verification failed!");
            }

            (setup_time, prove_time, verify_time)
        };

        println!("\n=== MPC Mode Stats ===");

        println!("Stats: {:#?}", MpcMultiNet::stats());

        MpcMultiNet::deinit();

        println!("\n=== MPC Mode Time Summary ===");
        println!("  n={}, m={}", n, m);
        println!("  {}: {:.4?}", method_name, permutation_time_elapsed);
        println!("  Proof setup:        {:.4?}", setup_time);
        println!("  Proof generation:   {:.4?}", prove_time);
        println!("  Proof verification: {:.4?}", verify_time);
    } else {
        let mut rng = test_rng();
        let (t_plain, f_plain) = prepare_input(n, m, &mut rng);

        let t: Vec<MF> = t_plain.iter().map(|&x| MF::Public(x)).collect();
        let f: Vec<MF> = f_plain.iter().map(|&x| MF::Public(x)).collect();

        let method_name = if opt.naive {
            "naive_secure_oblivious_lookup_permutation"
        } else {
            "secure_oblivious_lookup_permutation"
        };

        println!(
            "Running {} in local/public mode (n={}, m={})",
            method_name, n, m
        );

        let permutation_start_time = std::time::Instant::now();
        let permutation_timer =
            start_timer!(|| format!("{} local n={}, m={}", method_name, n, m));
        let (t_prime_mpc, f_prime_mpc) = if opt.naive {
            naive_secure_oblivious_lookup_permutation(n, m, &t, &f, opt.debug)
        } else {
            secure_oblivious_lookup_permutation(n, m, &t, &f, opt.debug)
        };
        end_timer!(permutation_timer);

        let permutation_time_elapsed = permutation_start_time.elapsed();

        let t_prime: Vec<F> = t_prime_mpc
            .iter()
            .map(|x| match x {
                MF::Public(val) => *val,
                _ => panic!("Expected Public values in local mode"),
            })
            .collect();
        let f_prime_unpadded: Vec<F> = f_prime_mpc
            .iter()
            .map(|x| match x {
                MF::Public(val) => *val,
                _ => panic!("Expected Public values in local mode"),
            })
            .collect();

        let f_prime_mpc = pad_f_prime_to_match_t_prime(f_prime_mpc, &t_prime_mpc, opt.debug);

        let f_prime: Vec<F> = f_prime_mpc
            .iter()
            .map(|x| match x {
                MF::Public(val) => *val,
                _ => panic!("Expected Public values in local mode"),
            })
            .collect();

        if opt.debug {
            println!(
                "t' length = {}, f' length = {} (padded)",
                t_prime.len(),
                f_prime.len()
            );
        }

        if opt.debug && n <= 16 {
            println!("\nDEBUG: After secure_oblivious_lookup_permutation:");
            println!("f (original queries) = {:?}", f_plain.iter().take(8).collect::<Vec<_>>());
            println!("f' (sorted queries, unpadded) = {:?}", f_prime_unpadded.iter().take(8).collect::<Vec<_>>());
            println!("f' (sorted queries, padded) = {:?}", f_prime.iter().take(8).collect::<Vec<_>>());
            println!("t (original table) = {:?}", t_plain.iter().take(8).collect::<Vec<_>>());
            println!("t' (permuted table) = {:?}", t_prime.iter().take(8).collect::<Vec<_>>());
        }

        if opt.debug {
            verify_lookup_permutation_outputs_in_plaintext(n, m, &t_plain, &f_plain, &t_prime, &f_prime_unpadded, opt.debug);
        }

        if opt.debug && n <= 16 {
            for (i, v) in t_prime.iter().enumerate() {
                println!("t'[{}] = {}", i, v);
            }
            for (i, v) in f_prime.iter().enumerate() {
                println!("f'[{}] = {}", i, v);
            }
        }

        let (setup_time, prove_time, verify_time) = {
            println!("\n=== Halo2 Lookup Argument PROOF Generation (Local) ===");
            
            let (config, ck, vk_pcs_local, setup_time) = setup_halo2_proof_system(n, m, opt.debug);
            
            let f_padded = pad_f_to_match_t(f.clone(), &t, &t_prime_mpc, opt.debug);
            let s_values = t.clone();
            let a_values = f_padded;
            
            let (revealed_proof, prove_time) = generate_halo2_proof(
                config.clone(),
                ck,
                a_values,
                s_values,
                t_prime_mpc.clone(),
                f_prime_mpc.clone(),
                opt.debug,
            );
            
            let (is_valid, verify_time) = verify_halo2_proof(config, vk_pcs_local, &revealed_proof, opt.debug);
            
            if !is_valid {
                eprintln!("WARNING: Proof verification failed!");
            }

            (setup_time, prove_time, verify_time)
        };

        println!("\n=== Local Mode Time Summary ===");        
        println!("  n={}, m={}", n, m);
        println!("  {}: {:.4?}", method_name, permutation_time_elapsed);
        println!("  Proof setup:        {:.4?}", setup_time);
        println!("  Proof generation:   {:.4?}", prove_time);
        println!("  Proof verification: {:.4?}", verify_time);
    }
}
