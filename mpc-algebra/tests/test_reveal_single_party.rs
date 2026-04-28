//! Test to demonstrate and verify the fix for reveal() in single-party mode
//! 
//! Before fix: Net::broadcast() returns an empty vector in single-party mode,
//! causing reveal() to sum an empty iterator, which returns zero.
//! 
//! After fix: When broadcast returns empty, reveal() returns the local value directly.

use ark_bls12_377::Fr;
use mpc_algebra::share::add::AdditiveFieldShare;
use mpc_algebra::Reveal;

#[test]
fn test_reveal_works_in_single_party_mode_after_fix() {
    println!("\n=== Test: reveal() works correctly in single-party mode after fix ===");
    
    // Create a test value
    let test_value = Fr::from(42u64);
    println!("Original value: {}", test_value);
    
    // Create a Shared value using from_add_shared
    let shared = AdditiveFieldShare::from_add_shared(test_value);
    println!("Shared value (stored locally): {}", shared.val);
    
    // Reveal the shared value
    // With the fix: when broadcast returns empty, return local value
    let revealed = shared.reveal();
    println!("Revealed value: {}", revealed);
    
    // After fix: revealed value should equal the original value
    assert_eq!(revealed, test_value, "After fix, reveal() returns the actual value in single-party mode");
    
    println!("✓ Test passes: reveal() now returns the correct value!");
    println!("  - Original value: {}", test_value);
    println!("  - Locally stored: {}", shared.val);
    println!("  - Revealed value: {} (correct!)", revealed);
    println!("\nThe fix:");
    println!("  if broadcast_results.is_empty() {{");
    println!("      self.val  // Return local value when no network");
    println!("  }} else {{");
    println!("      broadcast_results.into_iter().sum()");
    println!("  }}");
}

#[test]
fn test_reveal_works_for_public_after_fix() {
    println!("\n=== Test: reveal() works for Public values after fix ===");
    
    let test_value = Fr::from(100u64);
    
    // Create a Public value
    let public_share = AdditiveFieldShare::from_public(test_value);
    println!("Public share (king's value): {}", public_share.val);
    
    // Reveal should work for public values too
    let revealed = public_share.reveal();
    println!("Revealed value: {}", revealed);
    
    assert_eq!(revealed, test_value, "Public values also reveal correctly");
    println!("✓ Public values work correctly: {}", revealed);
}

#[test]
fn test_both_shared_and_public_work() {
    println!("\n=== Test: Both Shared and Public reveal correctly after fix ===");
    
    let test_value = Fr::from(999u64);
    
    // Test Shared
    let shared = AdditiveFieldShare::from_add_shared(test_value);
    let shared_revealed = shared.reveal();
    println!("Shared  - Original: {}, Revealed: {}", test_value, shared_revealed);
    assert_eq!(shared_revealed, test_value);
    
    // Test Public
    let public_share = AdditiveFieldShare::from_public(test_value);
    let public_revealed = public_share.reveal();
    println!("Public  - Original: {}, Revealed: {}", test_value, public_revealed);
    assert_eq!(public_revealed, test_value);
    
    println!("✓ Both Shared and Public reveal correctly!");
}
