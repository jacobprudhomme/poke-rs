use num_bigint::BigUint;

use crate::bn::BigNum;

// WARN: Assumes we are given an a priori invertible element as input
pub fn invert_element_mod(element: &BigNum, modulus: &BigNum) -> BigNum {
    let element = BigUint::from_bytes_le(element.as_le_bytes());
    let modulus = BigUint::from_bytes_le(modulus.as_le_bytes());

    let inverse = element.modinv(&modulus);
    let Some(inverse) = inverse else {
        unreachable!("We expect an invertible element as input");
    };

    BigNum::new(&inverse.to_u64_digits())
}
