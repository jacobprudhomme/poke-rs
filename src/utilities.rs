use num_bigint::BigUint;

use crate::bn::BigNum;

// WARN: Assumes we are given an a priori invertible element as input
pub fn invert_element_mod(element: &BigUint, modulus: &BigUint) -> BigNum {
    let inverse = element.modinv(modulus);
    let Some(inverse) = inverse else {
        unreachable!("We expect an invertible element as input");
    };

    BigNum::new(&inverse.to_u64_digits())
}
