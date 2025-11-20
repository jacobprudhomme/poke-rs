#![allow(non_upper_case_globals)]

use poke::{SUCCESS_RETVAL, solve_dlp_small_prime_power_order};
use rstest::rstest;

// p = 2 * 3^3 * 5^3 * 7^4 * 113 * 379 * 503 * 52837 + 1, such that |p| > 64 bits and Fp^2 has a multiplicative subgroup of order 5^3
const p: [u64; 2] = [0x0000000000000c3f, 0x0000000000000001];
fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = p,);

#[rstest]
fn test_simple_dlp() {
    // Generates the order-5^3 subgroup
    let generator = Fp2::const_decode_no_check(
        &[47, 47, 195, 76, 26, 17, 18, 144, 0],
        &[0; Fp::ENCODED_LENGTH],
    );
    // g^7
    let value = Fp2::const_decode_no_check(
        &[242, 200, 149, 124, 105, 197, 133, 202, 0],
        &[0; Fp::ENCODED_LENGTH],
    );

    assert_eq!(
        generator.pow(&[125], 7).equals(&Fp2::ONE),
        SUCCESS_RETVAL,
        "Generator doesn't have order 5^3"
    );

    let (x, ok) = solve_dlp_small_prime_power_order(&generator, &value, 5, 3);

    assert_eq!(ok, SUCCESS_RETVAL);
    assert_eq!(&x.repr, &[7]);
    assert_eq!(x.bitlen, 3);
}
#[rstest]
fn test_complex_dlp() {
    // Generates the order-7^4 subgroup
    let generator = Fp2::const_decode_no_check(
        &[82, 1, 230, 158, 56, 88, 5, 102, 0],
        &[0; Fp::ENCODED_LENGTH],
    );
    // g^501
    let value = Fp2::const_decode_no_check(
        &[148, 122, 116, 74, 152, 105, 50, 33, 0],
        &[0; Fp::ENCODED_LENGTH],
    );

    assert_eq!(
        generator.pow(&[97, 9], 12).equals(&Fp2::ONE),
        SUCCESS_RETVAL,
        "Generator doesn't have order 7^4"
    );

    let (x, ok) = solve_dlp_small_prime_power_order(&generator, &value, 7, 4);

    assert_eq!(ok, SUCCESS_RETVAL);
    assert_eq!(&x.repr, &[255, 1]);
    assert_eq!(x.bitlen, 9);
}
