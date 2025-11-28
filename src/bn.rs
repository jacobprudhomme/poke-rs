use isogeny::utilities::bn::{
    bn_bit_length_vartime, factorisation_to_bn_vartime, prime_power_to_bn_vartime,
};

#[derive(Debug, PartialEq)]
pub struct BigNum {
    pub repr: Vec<u8>,
    pub bitlen: usize,
}

impl BigNum {
    pub fn new(le_words: &[u64]) -> Self {
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

        Self {
            repr: le_bytes,
            bitlen,
        }
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
        self.repr
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

    pub fn bits(&self) -> usize {
        self.bitlen
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
