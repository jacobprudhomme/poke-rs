use num_bigint::BigUint;

use crate::bn::BigNum;

// WARN: Assumes we are given an a priori invertible element as input
pub fn invert_element_mod(element: &BigUint, modulus: &BigUint) -> BigNum {
    let inverse = element.modinv(modulus);
    let Some(inverse) = inverse else {
        unreachable!("We expect an invertible element as input");
    };

    let bitlen = inverse
        .bits()
        .try_into()
        .expect("Size in bits of the inverse scalar is too big to fit in a usize (we do not ever expect this to happen)");
    let repr = inverse.to_bytes_le();

    BigNum { repr, bitlen }
}
