use isogeny::utilities::bn::bn_bit_length_vartime;

#[derive(Debug, PartialEq)]
pub struct BigNum {
    pub repr: Vec<u8>,
    pub bitlen: usize,
}

pub fn byte_bn_to_word_bn(bytes: &BigNum) -> Vec<u64> {
    bytes
        .repr
        .chunks(8)
        .map(|word_bytes| {
            let mut word_bytes = word_bytes.to_vec();
            word_bytes.resize(8, 0);
            u64::from_le_bytes(word_bytes.try_into().unwrap())
        })
        .collect()
}

pub fn word_bn_to_byte_bn(words: &[u64]) -> BigNum {
    let bitlen = bn_bit_length_vartime(words);
    let mut bytes = words
        .iter()
        .flat_map(|word| word.to_le_bytes())
        .collect::<Vec<_>>();
    while bytes.len() > 1
        && let Some(&last_byte) = bytes.last()
        && last_byte == 0
    {
        bytes.pop();
    }

    BigNum {
        repr: bytes,
        bitlen,
    }
}

#[cfg(test)]
mod tests {
    use rand::RngCore as _;
    use rstest::rstest;

    use super::{BigNum, byte_bn_to_word_bn, word_bn_to_byte_bn};

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

        let bn_words = byte_bn_to_word_bn(&bn);
        let bn_bytes_roundtrip = word_bn_to_byte_bn(&bn_words);

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

        let bn = word_bn_to_byte_bn(&bn_words);
        let bn_words_roundtrip = byte_bn_to_word_bn(&bn);

        assert_eq!(bn_words, bn_words_roundtrip);
    }
}
