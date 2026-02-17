use core::marker::PhantomData;

use isogeny::utilities::bn::bn_mul_vartime;
use rug::{Integer, integer::Order};

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
    pub fn invert_mod_vartime<const NUM_WORDS_MOD: usize>(
        &self,
        modulus: &BigNum<NUM_WORDS_MOD>,
    ) -> BigNum<NUM_WORDS_MOD> {
        let this = Integer::from_digits(self.as_le_words(), Order::Lsf);
        let modulus = Integer::from_digits(modulus.as_le_words(), Order::Lsf);

        let Ok(inverse) = this.invert(&modulus) else {
            unreachable!("Input should be a unit with respect to the modulus");
        };

        BigNum::new(inverse.as_limbs())
    }
}

pub fn crt2<
    const NUM_WORDS_X: usize,
    const NUM_WORDS_Y: usize,
    const NUM_WORDS_XY: usize,
    const NUM_WORDS_XXY: usize,
    const NUM_WORDS_XYY: usize,
>(
    residues: (&BigNum<NUM_WORDS_X>, &BigNum<NUM_WORDS_Y>),
    partial_moduli: (&BigNum<NUM_WORDS_X>, &BigNum<NUM_WORDS_Y>),
    inv_partial_moduli: (&BigNum<NUM_WORDS_X>, &BigNum<NUM_WORDS_Y>),
    full_modulus: &BigNum<NUM_WORDS_XY>,
    _intermediate_bignum_sizes: PhantomData<([(); NUM_WORDS_XXY], [(); NUM_WORDS_XYY])>,
) -> BigNum<NUM_WORDS_XY> {
    let mut result = BigNum::zero();

    // FIXME: Is it possible that the result goes beyond its size bounds when adding the reduced partial results?
    // Or does BigNUM<NUM_WORDS_XY> have enough room for the result to grow a bit, for all concrete instantiations?

    // Solution mod x^ex
    let partial_result: BigNum<NUM_WORDS_XXY> =
        residues.0.widen() * inv_partial_moduli.0.widen() * partial_moduli.1.widen();
    result += partial_result.reduce_mod(full_modulus);

    // Solution mod y^ey
    let partial_result: BigNum<NUM_WORDS_XYY> =
        residues.1.widen() * inv_partial_moduli.1.widen() * partial_moduli.0.widen();
    result = (result + partial_result.reduce_mod(full_modulus)).reduce_mod(full_modulus);

    result
}

pub fn crt3<
    const NUM_WORDS_X: usize,
    const NUM_WORDS_Y: usize,
    const NUM_WORDS_Z: usize,
    const NUM_WORDS_XY: usize,
    const NUM_WORDS_XZ: usize,
    const NUM_WORDS_YZ: usize,
    const NUM_WORDS_XYZ: usize,
    const NUM_WORDS_XXYZ: usize,
    const NUM_WORDS_XYYZ: usize,
    const NUM_WORDS_XYZZ: usize,
>(
    residues: (
        &BigNum<NUM_WORDS_X>,
        &BigNum<NUM_WORDS_Y>,
        &BigNum<NUM_WORDS_Z>,
    ),
    duals_of_partial_moduli: (
        &BigNum<NUM_WORDS_YZ>,
        &BigNum<NUM_WORDS_XZ>,
        &BigNum<NUM_WORDS_XY>,
    ),
    inv_duals_of_partial_moduli: (
        &BigNum<NUM_WORDS_X>,
        &BigNum<NUM_WORDS_Y>,
        &BigNum<NUM_WORDS_Z>,
    ),
    full_modulus: &BigNum<NUM_WORDS_XYZ>,
    _intermediate_bignum_sizes: PhantomData<(
        [(); NUM_WORDS_XXYZ],
        [(); NUM_WORDS_XYYZ],
        [(); NUM_WORDS_XYZZ],
    )>,
) -> BigNum<NUM_WORDS_XYZ> {
    let mut result = BigNum::zero();

    // FIXME: Is it possible that the result goes beyond its size bounds when adding the reduced partial results?
    // Or does BigNUM<NUM_WORDS_XYZ> have enough room for the result to grow a bit, for all concrete instantiations?

    // Solution mod x^ex
    let partial_result: BigNum<NUM_WORDS_XXYZ> = residues.0.widen()
        * inv_duals_of_partial_moduli.0.widen()
        * duals_of_partial_moduli.0.widen();
    result += partial_result.reduce_mod(full_modulus);

    // Solution mod y^ey
    let partial_result: BigNum<NUM_WORDS_XYYZ> = residues.1.widen()
        * inv_duals_of_partial_moduli.1.widen()
        * duals_of_partial_moduli.1.widen();
    result = (result + partial_result.reduce_mod(full_modulus)).reduce_mod(full_modulus);

    // Solution mod z^ez
    let partial_result: BigNum<NUM_WORDS_XYZZ> = residues.2.widen()
        * inv_duals_of_partial_moduli.2.widen()
        * duals_of_partial_moduli.2.widen();
    result = (result + partial_result.reduce_mod(full_modulus)).reduce_mod(full_modulus);

    result
}
