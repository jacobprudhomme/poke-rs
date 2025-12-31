use core::{
    array, cmp, fmt,
    ops::{Add, AddAssign, Mul, MulAssign},
};
use std::ops::Sub;

use isogeny::utilities::bn::{
    bn_add_vartime, bn_bit_length_vartime, bn_mul_by_u64_vartime, bn_mul_vartime,
    bn_sub_into_vartime, factorisation_to_bn_vartime, prime_power_to_bn_vartime,
};
use num_bigint::BigUint;

#[derive(Debug, PartialEq)]
pub struct BigNumArb {
    repr: Vec<u64>,
    bitlen: usize,
}

#[derive(Debug)]
pub struct BigNum<const NUM_WORDS: usize> {
    repr: [u64; NUM_WORDS],
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

impl BigNumArb {
    pub fn new(words: &[u64]) -> Self {
        let bitlen = bn_bit_length_vartime(&words);

        Self {
            repr: words.to_owned(),
            bitlen,
        }
    }

    pub fn from_le_bytes(bytes: &[u8]) -> Self {
        Self::new(&le_bytes_to_le_words(bytes))
    }

    pub fn as_le_words<'a>(&'a self) -> &'a [u64] {
        &self.repr
    }
}

impl<const NUM_WORDS: usize> BigNum<NUM_WORDS> {
    pub fn zero() -> Self {
        Self {
            repr: [0; NUM_WORDS],
            bitlen: 0,
        }
    }

    pub fn one() -> Self {
        let mut words = [1; NUM_WORDS];
        words[0] = 1;

        Self {
            repr: words,
            bitlen: 1,
        }
    }

    // WARN: 0-extends if le_words.len() < NUM_WORDS, truncates if le_words.len() > NUM_WORDS
    pub fn new(words: &[u64]) -> Self {
        let least_num_of_words = cmp::min(words.len(), NUM_WORDS);

        let mut repr = [0u64; NUM_WORDS];
        // SAFETY: Source and destination will necessarily be the same size, and the
        // length of the copied slice is less than or equal to the length of both
        repr[..least_num_of_words].copy_from_slice(&words[..least_num_of_words]);

        let bitlen = bn_bit_length_vartime(&repr);

        Self { repr, bitlen }
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

    pub fn to_le_bytes(&self) -> Vec<u8> {
        le_words_to_le_bytes(&self.repr).0
    }

    pub fn as_le_words<'a>(&'a self) -> &'a [u64] {
        &self.repr
    }

    pub fn nbits(&self) -> usize {
        self.bitlen
    }

    pub fn widening_mul<const NUM_WORDS_RHS: usize>(
        &self,
        rhs: &BigNum<NUM_WORDS_RHS>,
    ) -> BigNum<{ NUM_WORDS + NUM_WORDS_RHS }> {
        BigNum::new(&bn_mul_vartime(&self.repr, &rhs.repr))
    }

    pub fn truncate<const NEW_NUM_WORDS: usize>(&self) -> BigNum<NEW_NUM_WORDS> {
        let new_max_bitlen = NEW_NUM_WORDS * 64;
        let bitlen = cmp::min(self.bitlen, new_max_bitlen);

        BigNum {
            repr: array::from_fn(|i| self.repr[i]),
            bitlen,
        }
    }
}

// WARN: Since addition is expected to have size max(lhs.len(), rhs.len()) + 1,
// and multiplication is expected to have size lhs.len() + rhs.len(), all of the
// following impls and can easily truncate the result. Use for results that are
// expected to be within the appropriate word size.

impl<const NUM_WORDS: usize> Add<&BigNum<NUM_WORDS>> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn add(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Add<BigNum<NUM_WORDS>> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn add(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Add<&BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn add(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Add<BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn add(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_add_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Sub<&BigNum<NUM_WORDS>> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn sub(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        let mut words = self.repr;
        bn_sub_into_vartime(&mut words, &rhs.repr);
        BigNum::new(&words)
    }
}

impl<const NUM_WORDS: usize> Sub<BigNum<NUM_WORDS>> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn sub(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        let mut words = self.repr;
        bn_sub_into_vartime(&mut words, &rhs.repr);
        BigNum::new(&words)
    }
}

impl<const NUM_WORDS: usize> Sub<&BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn sub(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        let mut words = self.repr;
        bn_sub_into_vartime(&mut words, &rhs.repr);
        BigNum::new(&words)
    }
}

impl<const NUM_WORDS: usize> Sub<BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn sub(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        let mut words = self.repr;
        bn_sub_into_vartime(&mut words, &rhs.repr);
        BigNum::new(&words)
    }
}

impl<const NUM_WORDS: usize> Mul<&BigNum<NUM_WORDS>> for u64 {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&rhs.repr, self))
    }
}

impl<const NUM_WORDS: usize> Mul<BigNum<NUM_WORDS>> for u64 {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&rhs.repr, self))
    }
}

impl<const NUM_WORDS: usize> Mul<u64> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: u64) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&self.repr, rhs))
    }
}

impl<const NUM_WORDS: usize> Mul<u64> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: u64) -> Self::Output {
        BigNum::new(&bn_mul_by_u64_vartime(&self.repr, rhs))
    }
}

impl<const NUM_WORDS: usize> Mul<&BigNum<NUM_WORDS>> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Mul<BigNum<NUM_WORDS>> for &BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Mul<&BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: &BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> Mul<BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    type Output = BigNum<NUM_WORDS>;

    fn mul(self, rhs: BigNum<NUM_WORDS>) -> Self::Output {
        BigNum::new(&bn_mul_vartime(&self.repr, &rhs.repr))
    }
}

impl<const NUM_WORDS: usize> AddAssign<&BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    fn add_assign(&mut self, rhs: &BigNum<NUM_WORDS>) {
        let words = bn_add_vartime(&self.repr, &rhs.repr);
        let bn = BigNum::new(&words);
        self.repr = bn.repr;
        self.bitlen = bn.bitlen;
    }
}

impl<const NUM_WORDS: usize> AddAssign<BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    fn add_assign(&mut self, rhs: BigNum<NUM_WORDS>) {
        let words = bn_add_vartime(&self.repr, &rhs.repr);
        let bn = BigNum::new(&words);
        self.repr = bn.repr;
        self.bitlen = bn.bitlen;
    }
}

impl<const NUM_WORDS: usize> MulAssign<&BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    fn mul_assign(&mut self, rhs: &BigNum<NUM_WORDS>) {
        let words = bn_mul_vartime(&self.repr, &rhs.repr);
        let bn = BigNum::new(&words);
        self.repr = bn.repr;
        self.bitlen = bn.bitlen;
    }
}

impl<const NUM_WORDS: usize> MulAssign<BigNum<NUM_WORDS>> for BigNum<NUM_WORDS> {
    fn mul_assign(&mut self, rhs: BigNum<NUM_WORDS>) {
        let words = bn_mul_vartime(&self.repr, &rhs.repr);
        let bn = BigNum::new(&words);
        self.repr = bn.repr;
        self.bitlen = bn.bitlen;
    }
}

impl<const NUM_WORDS: usize> fmt::Display for BigNum<NUM_WORDS> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        BigUint::from_bytes_le(&self.to_le_bytes()).fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use rand::RngCore as _;
    use rstest::rstest;

    use super::BigNumArb;

    #[rstest]
    fn word_bn_is_inverse_of_byte_bn() {
        let mut rng = rand::rng();
        let mut bn_bytes = vec![0u8; 43];
        rng.fill_bytes(&mut bn_bytes[..42]);
        bn_bytes[42] = 1;

        let bn = BigNumArb::from_le_bytes(&bn_bytes);
        let bn_words = bn.as_le_words();
        let bn_bytes_roundtrip = BigNumArb::new(bn_words);

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
        let bn = BigNumArb::new(&bn_words);
        let bn_words_roundtrip = bn.as_le_words();

        assert_eq!(bn_words, bn_words_roundtrip);
    }
}
