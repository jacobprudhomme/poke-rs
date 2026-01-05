use isogeny::utilities::bn::bn_mul_vartime;
use num_bigint::BigUint;

use crate::bn::BigNum;

impl<const NUM_WORDS: usize> BigNum<NUM_WORDS> {
    pub fn mul_mod_power_of_two<const NUM_WORDS_RHS: usize, const NUM_WORDS_MOD: usize>(
        &self,
        rhs: &BigNum<NUM_WORDS_RHS>,
        modulus: &BigNum<NUM_WORDS_MOD>,
    ) -> Self {
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

    pub fn reduce_mod<const NUM_WORDS_MOD: usize>(
        &self,
        modulus: &BigNum<NUM_WORDS_MOD>,
    ) -> BigNum<NUM_WORDS_MOD> {
        let this = BigUint::from_bytes_le(&self.to_le_bytes());
        let modulus = BigUint::from_bytes_le(&modulus.to_le_bytes());

        let result = this % modulus;
        BigNum::new(&result.to_u64_digits())
    }

    // WARN: Assumes n is a unit with respect to the modulus
    pub fn invert_mod<const NUM_WORDS_MOD: usize>(
        &self,
        modulus: &BigNum<NUM_WORDS_MOD>,
    ) -> BigNum<NUM_WORDS_MOD> {
        let this = BigUint::from_bytes_le(&self.to_le_bytes());
        let modulus = BigUint::from_bytes_le(&modulus.to_le_bytes());

        let Some(inverse) = this.modinv(&modulus) else {
            unreachable!("Input should have been a unit with respect to the modulus");
        };

        BigNum::new(&inverse.to_u64_digits())
    }
}

pub fn crt_mod_powers_of_two_three_five<
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_235: usize,
>(
    residues: (
        &BigNum<NUM_WORDS_2>,
        &BigNum<NUM_WORDS_3>,
        &BigNum<NUM_WORDS_5>,
    ),
    partial_moduli: (
        &BigNum<NUM_WORDS_2>,
        &BigNum<NUM_WORDS_3>,
        &BigNum<NUM_WORDS_5>,
    ),
    partial_modulus_products: (
        &BigNum<NUM_WORDS_35>,
        &BigNum<NUM_WORDS_25>,
        &BigNum<NUM_WORDS_23>,
    ),
    full_modulus: &BigNum<NUM_WORDS_235>,
) -> BigNum<NUM_WORDS_235> {
    let mut result = BigNum::zero();

    // FIXME: Is it possible that any of these go beyond the prescribed fixed word-size?
    result += residues.0.widen()
        * partial_modulus_products
            .0
            .invert_mod(partial_moduli.0)
            .widen()
        * partial_modulus_products.0.widen();
    result = result.reduce_mod(full_modulus);
    result += residues.1.widen()
        * partial_modulus_products
            .1
            .invert_mod(partial_moduli.1)
            .widen()
        * partial_modulus_products.1.widen();
    result = result.reduce_mod(full_modulus);
    result += residues.2.widen()
        * partial_modulus_products
            .2
            .invert_mod(partial_moduli.2)
            .widen()
        * partial_modulus_products.2.widen();
    result = result.reduce_mod(full_modulus);

    result
}
