use isogeny::utilities::bn::bn_mul_vartime;

use crate::bn::BigNum;

impl<const NUM_WORDS: usize> BigNum<NUM_WORDS> {
    pub fn mul_mod_power_of_two(&self, rhs: &Self, modulus: &Self) -> Self {
        let mut result = bn_mul_vartime(&self.as_le_words(), &rhs.as_le_words());

        let mut ctl = u64::MAX;
        let mut mask = [0u64; NUM_WORDS];
        for (i, modulus_word) in modulus.as_le_words().iter().enumerate() {
            mask[i] = ctl ^ (!modulus_word).wrapping_add(1);

            // 0 if modulus_word is 0, 1 otherwise
            let modulus_word_is_nonzero =
                (modulus_word | modulus_word.wrapping_neg()) >> (u64::BITS - 1);
            // Expand bit to entire 64-bit word
            let modulus_word_is_nonzero = modulus_word_is_nonzero.wrapping_neg();
            ctl ^= modulus_word_is_nonzero;
        }

        for (result_word, mask) in result.iter_mut().zip(mask) {
            *result_word &= mask;
        }

        BigNum::new(&result)
    }
}
