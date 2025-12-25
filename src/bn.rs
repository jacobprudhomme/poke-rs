use core::{
    fmt,
    ops::{Add, Mul},
};
use std::ops::{AddAssign, MulAssign};

use isogeny::utilities::bn::{
    bn_add_vartime, bn_bit_length_vartime, bn_from_le_bytes, bn_mul_by_u64_vartime, bn_mul_vartime,
    factorisation_to_bn_vartime, prime_power_to_bn_vartime,
};
use num_bigint::BigUint;

#[derive(Debug, PartialEq)]
pub struct BigNum {
    repr: Vec<u8>,
    bitlen: usize,
}

fn le_bytes_to_le_words(bytes: &[u8]) -> Vec<u64> {
    bytes
        .chunks(8)
        .map(|word_bytes| {
            let mut word_bytes = word_bytes.to_vec();
            word_bytes.resize(8, 0);
            let Ok(word_bytes) = word_bytes.try_into() else {
                unreachable!("Need 8 bytes to form a u64 (we never expect to reach here)")
            };
            u64::from_le_bytes(word_bytes)
        })
        .collect()
}

fn le_words_to_le_bytes(le_words: &[u64]) -> (Vec<u8>, usize) {
    let mut le_bytes = le_words
        .iter()
        .flat_map(|word| word.to_le_bytes())
        .collect::<Vec<_>>();
    while le_bytes.len() > 1
        && let Some(&last_byte) = le_bytes.last()
        && last_byte == 0
    {
        le_bytes.pop();
    }
    let bitlen = bn_bit_length_vartime(le_words);

    (le_bytes, bitlen)
}

// WARN: all of these functions are vartime
impl BigNum {
    pub fn zero() -> Self {
        Self {
            repr: vec![0],
            bitlen: 0,
        }
    }

    pub fn one() -> Self {
        Self {
            repr: vec![1],
            bitlen: 1,
        }
    }

    pub fn new(le_words: &[u64]) -> Self {
        let (le_bytes, bitlen) = le_words_to_le_bytes(le_words);

        Self {
            repr: le_bytes,
            bitlen,
        }
    }

    pub fn from_le_bytes(bytes: &[u8]) -> Self {
        Self::new(&le_bytes_to_le_words(bytes))
    }

    pub fn from_prime(p: usize) -> Self {
        Self::new(&prime_power_to_bn_vartime(p, 1))
    }

    pub fn from_prime_power(p: usize, exp: usize) -> Self {
        Self::new(&prime_power_to_bn_vartime(p, exp))
    }

    pub fn from_prime_factors(prime_factorization: &[(usize, usize)]) -> Self {
        Self::new(&factorisation_to_bn_vartime(prime_factorization))
    }

    pub fn as_le_bytes<'a>(&'a self) -> &'a [u8] {
        &self.repr
    }

    pub fn to_le_words(&self) -> Vec<u64> {
        bn_from_le_bytes(&self.repr, self.bitlen)
    }

    pub fn nbits(&self) -> usize {
        self.bitlen
    }
}

impl Add<&BigNum> for &BigNum {
    type Output = BigNum;

    fn add(self, rhs: &BigNum) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl Add<BigNum> for &BigNum {
    type Output = BigNum;

    fn add(self, rhs: BigNum) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl Add<&BigNum> for BigNum {
    type Output = BigNum;

    fn add(self, rhs: &BigNum) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl Add<BigNum> for BigNum {
    type Output = BigNum;

    fn add(self, rhs: BigNum) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl AddAssign<&BigNum> for BigNum {
    fn add_assign(&mut self, rhs: &BigNum) {
        let le_words = bn_add_vartime(&self.to_le_words(), &rhs.to_le_words());
        (self.repr, self.bitlen) = le_words_to_le_bytes(&le_words);
    }
}

impl AddAssign<BigNum> for BigNum {
    fn add_assign(&mut self, rhs: BigNum) {
        let le_words = bn_add_vartime(&self.to_le_words(), &rhs.to_le_words());
        (self.repr, self.bitlen) = le_words_to_le_bytes(&le_words);
    }
}

impl Mul<&BigNum> for u64 {
    type Output = BigNum;

    fn mul(self, rhs: &BigNum) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&rhs.to_le_words(), self))
    }
}

impl Mul<BigNum> for u64 {
    type Output = BigNum;

    fn mul(self, rhs: BigNum) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&rhs.to_le_words(), self))
    }
}

impl Mul<u64> for &BigNum {
    type Output = BigNum;

    fn mul(self, rhs: u64) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&self.to_le_words(), rhs))
    }
}

impl Mul<u64> for BigNum {
    type Output = BigNum;

    fn mul(self, rhs: u64) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&self.to_le_words(), rhs))
    }
}

impl Mul<&BigNum> for &BigNum {
    type Output = BigNum;

    fn mul(self, rhs: &BigNum) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl Mul<BigNum> for &BigNum {
    type Output = BigNum;

    fn mul(self, rhs: BigNum) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl Mul<&BigNum> for BigNum {
    type Output = BigNum;

    fn mul(self, rhs: &BigNum) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl Mul<BigNum> for BigNum {
    type Output = BigNum;

    fn mul(self, rhs: BigNum) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.to_le_words(), &rhs.to_le_words()))
    }
}

impl MulAssign<&BigNum> for BigNum {
    fn mul_assign(&mut self, rhs: &BigNum) {
        let le_words = bn_mul_vartime(&self.to_le_words(), &rhs.to_le_words());
        (self.repr, self.bitlen) = le_words_to_le_bytes(&le_words);
    }
}

impl MulAssign<BigNum> for BigNum {
    fn mul_assign(&mut self, rhs: BigNum) {
        let le_words = bn_mul_vartime(&self.to_le_words(), &rhs.to_le_words());
        (self.repr, self.bitlen) = le_words_to_le_bytes(&le_words);
    }
}

impl fmt::Display for BigNum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        BigUint::from_bytes_le(self.as_le_bytes()).fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use rand::RngCore as _;
    use rstest::rstest;

    use super::BigNum;

    #[rstest]
    fn word_bn_is_inverse_of_byte_bn() {
        let mut rng = rand::rng();
        let mut bn_bytes = vec![0u8; 43];
        rng.fill_bytes(&mut bn_bytes[..42]);
        bn_bytes[42] = 1;
        let bn = BigNum {
            repr: bn_bytes,
            bitlen: 42 * 8 + 1,
        };

        let bn_words = bn.to_le_words();
        let bn_bytes_roundtrip = BigNum::new(&bn_words);

        assert_eq!(bn, bn_bytes_roundtrip);
    }

    #[rstest]
    fn byte_bn_is_inverse_of_word_bn() {
        let mut rng = rand::rng();
        let mut bn_bytes = [0u8; 48];
        rng.fill_bytes(&mut bn_bytes);
        for i in 43..48 {
            bn_bytes[i] = 0;
        }
        let (_, bn_words, _) = unsafe { bn_bytes.align_to::<u64>() };

        let bn = BigNum::new(&bn_words);
        let bn_words_roundtrip = bn.to_le_words();

        assert_eq!(bn_words, bn_words_roundtrip);
    }
}
