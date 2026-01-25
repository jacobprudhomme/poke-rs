use core::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;

use crate::{
    FAILURE_RETVAL, SUCCESS_RETVAL,
    bn::BigNum,
    modular::{crt2, crt3},
};

// Works for primes < 256. Assumes generator `generator` generates the subgroup.
pub fn solve_dlp_small_prime_order<Fp2: Fp2Trait>(
    generator: &Fp2,
    value: &Fp2,
    p: u8,
) -> (u8, u32) {
    let mut retval = FAILURE_RETVAL;

    let mut element = Fp2::ONE;
    let mut result = 0;
    for i in 0..p {
        // Only set this if we have not already found the log
        let found_log = !retval & value.equals(&element);
        result |= i & (found_log as u8);
        retval |= found_log;

        element *= *generator;
    }

    (result, retval)
}

pub fn solve_dlp_small_prime_power_order<Fp2: Fp2Trait, const NUM_WORDS_ORDER: usize>(
    generator: &Fp2,
    value: &Fp2,
    p: u8,
    e: usize,
    p_adic_basis: &[BigNum<NUM_WORDS_ORDER>],
) -> (BigNum<NUM_WORDS_ORDER>, u32) {
    let mut retval = SUCCESS_RETVAL;

    let prime_order_subgroup_generator = generator.pow(
        &p_adic_basis[e - 1].to_le_bytes(),
        p_adic_basis[e - 1].nbits(),
    );

    let mut partial_sum = BigNum::zero();
    for i in 0..e {
        let r = *value
            * generator
                .pow(&partial_sum.to_le_bytes(), partial_sum.nbits())
                .invert(); // FIXME: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            &p_adic_basis[e - i - 1].to_le_bytes(),
            p_adic_basis[e - i - 1].nbits(),
        );

        let (x, ok) = solve_dlp_small_prime_order(&prime_order_subgroup_generator, &u, p);
        partial_sum += (x as u64) * &p_adic_basis[i];
        retval &= ok;
    }

    (partial_sum, retval)
}

pub fn solve_dlp_order_powers_of_two_three<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_233: usize,
>(
    generator: &Fp2,
    value: &Fp2,
    es: (usize, usize),
    prime_power_orders: (&BigNum<NUM_WORDS_2>, &BigNum<NUM_WORDS_3>),
    inv_prime_power_orders: (&BigNum<NUM_WORDS_2>, &BigNum<NUM_WORDS_3>),
    full_product_of_prime_power_orders: &BigNum<NUM_WORDS_23>,
    p_adic_bases: (&[BigNum<NUM_WORDS_2>], &[BigNum<NUM_WORDS_3>]),
    intermediate_bignum_sizes: PhantomData<([(); NUM_WORDS_223], [(); NUM_WORDS_233])>,
) -> (BigNum<NUM_WORDS_23>, u32) {
    let mut retval = SUCCESS_RETVAL;

    let generator_of_power_of_two_subgroup = generator.pow(
        &prime_power_orders.1.to_le_bytes(),
        prime_power_orders.1.nbits(),
    );
    let value_in_power_of_two_subgroup = value.pow(
        &prime_power_orders.1.to_le_bytes(),
        prime_power_orders.1.nbits(),
    );
    let (result_mod_power_of_two, ok) = solve_dlp_small_prime_power_order(
        &generator_of_power_of_two_subgroup,
        &value_in_power_of_two_subgroup,
        2,
        es.0,
        p_adic_bases.0,
    );
    retval &= ok;

    let generator_of_power_of_three_subgroup = generator.pow(
        &prime_power_orders.0.to_le_bytes(),
        prime_power_orders.0.nbits(),
    );
    let value_in_power_of_three_subgroup = value.pow(
        &prime_power_orders.0.to_le_bytes(),
        prime_power_orders.0.nbits(),
    );
    let (result_mod_power_of_three, ok) = solve_dlp_small_prime_power_order(
        &generator_of_power_of_three_subgroup,
        &value_in_power_of_three_subgroup,
        3,
        es.1,
        p_adic_bases.1,
    );
    retval &= ok;

    let result = crt2(
        (&result_mod_power_of_two, &result_mod_power_of_three),
        prime_power_orders,
        inv_prime_power_orders,
        full_product_of_prime_power_orders,
        intermediate_bignum_sizes,
    );

    (result, retval)
}

pub fn solve_dlp_order_powers_of_two_three_five<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_235: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
>(
    generator: &Fp2,
    value: &Fp2,
    es: (usize, usize, usize),
    partial_products_of_prime_power_orders: (
        &BigNum<NUM_WORDS_35>,
        &BigNum<NUM_WORDS_25>,
        &BigNum<NUM_WORDS_23>,
    ),
    inv_partial_products_of_prime_power_orders: (
        &BigNum<NUM_WORDS_2>,
        &BigNum<NUM_WORDS_3>,
        &BigNum<NUM_WORDS_5>,
    ),
    full_product_of_prime_power_orders: &BigNum<NUM_WORDS_235>,
    p_adic_bases: (
        &[BigNum<NUM_WORDS_2>],
        &[BigNum<NUM_WORDS_3>],
        &[BigNum<NUM_WORDS_5>],
    ),
    intermediate_bignum_sizes: PhantomData<(
        [(); NUM_WORDS_2235],
        [(); NUM_WORDS_2335],
        [(); NUM_WORDS_2355],
    )>,
) -> (BigNum<NUM_WORDS_235>, u32) {
    let mut retval = SUCCESS_RETVAL;

    let generator_of_power_of_two_subgroup = generator.pow(
        &partial_products_of_prime_power_orders.0.to_le_bytes(),
        partial_products_of_prime_power_orders.0.nbits(),
    );
    let value_in_power_of_two_subgroup = value.pow(
        &partial_products_of_prime_power_orders.0.to_le_bytes(),
        partial_products_of_prime_power_orders.0.nbits(),
    );
    let (result_mod_power_of_two, ok) = solve_dlp_small_prime_power_order(
        &generator_of_power_of_two_subgroup,
        &value_in_power_of_two_subgroup,
        2,
        es.0,
        p_adic_bases.0,
    );
    retval &= ok;

    let generator_of_power_of_three_subgroup = generator.pow(
        &partial_products_of_prime_power_orders.1.to_le_bytes(),
        partial_products_of_prime_power_orders.1.nbits(),
    );
    let value_in_power_of_three_subgroup = value.pow(
        &partial_products_of_prime_power_orders.1.to_le_bytes(),
        partial_products_of_prime_power_orders.1.nbits(),
    );
    let (result_mod_power_of_three, ok) = solve_dlp_small_prime_power_order(
        &generator_of_power_of_three_subgroup,
        &value_in_power_of_three_subgroup,
        3,
        es.1,
        p_adic_bases.1,
    );
    retval &= ok;

    let generator_of_power_of_five_subgroup = generator.pow(
        &partial_products_of_prime_power_orders.2.to_le_bytes(),
        partial_products_of_prime_power_orders.2.nbits(),
    );
    let value_in_power_of_five_subgroup = value.pow(
        &partial_products_of_prime_power_orders.2.to_le_bytes(),
        partial_products_of_prime_power_orders.2.nbits(),
    );
    let (result_mod_power_of_five, ok) = solve_dlp_small_prime_power_order(
        &generator_of_power_of_five_subgroup,
        &value_in_power_of_five_subgroup,
        5,
        es.2,
        p_adic_bases.2,
    );
    retval &= ok;

    let result = crt3(
        (
            &result_mod_power_of_two,
            &result_mod_power_of_three,
            &result_mod_power_of_five,
        ),
        partial_products_of_prime_power_orders,
        inv_partial_products_of_prime_power_orders,
        full_product_of_prime_power_orders,
        intermediate_bignum_sizes,
    );

    (result, retval)
}

pub fn solve_dlp_order_five<Fp2: Fp2Trait>(generator: &Fp2, value: &Fp2) -> (u8, u32) {
    let mut retval = FAILURE_RETVAL;

    let mut element = Fp2::ONE;
    let mut result = 0;

    let found_log = value.equals(&element);
    result |= 0 & (found_log as u8);
    retval |= found_log;

    element *= *generator;

    let found_log = value.equals(&element);
    result |= 1 & (found_log as u8);
    retval |= found_log;

    let mut element2 = element * *generator;

    let found_log = value.equals(&element2);
    result |= 2 & (found_log as u8);
    retval |= found_log;

    element2.set_conjugate();

    let found_log = value.equals(&element2);
    result |= 3 & (found_log as u8);
    retval |= found_log;

    element.set_conjugate();

    let found_log = value.equals(&element);
    result |= 4 & (found_log as u8);
    retval |= found_log;

    (result, retval)
}

pub fn solve_dlp_order_power_of_five<Fp2: Fp2Trait, const NUM_WORDS_ORDER: usize>(
    generator: &Fp2,
    value: &Fp2,
    e: usize,
    five_adic_basis: &[BigNum<NUM_WORDS_ORDER>],
) -> (BigNum<NUM_WORDS_ORDER>, u32) {
    let mut retval = SUCCESS_RETVAL;

    let prime_order_subgroup_generator = generator.pow(
        &five_adic_basis[e - 1].to_le_bytes(),
        five_adic_basis[e - 1].nbits(),
    );
    assert_eq!(
        prime_order_subgroup_generator.equals(&Fp2::ONE),
        FAILURE_RETVAL,
        "g has order < 5^{}",
        e,
    );

    let mut partial_sum = BigNum::zero();
    for i in 0..e {
        let r = *value
            * generator
                .pow(&partial_sum.to_le_bytes(), partial_sum.nbits())
                .invert(); // FIXME: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            &five_adic_basis[e - i - 1].to_le_bytes(),
            five_adic_basis[e - i - 1].nbits(),
        );

        let (x, ok) = solve_dlp_order_five(&prime_order_subgroup_generator, &u);
        partial_sum += (x as u64) * &five_adic_basis[i];
        retval &= ok;
    }

    (partial_sum, retval)
}

pub fn solve_dlp_order_five_explicit_subgroup<Fp2: Fp2Trait>(
    subgroup: &[Fp2; 5],
    value: &Fp2,
) -> (u8, u32) {
    let mut retval = FAILURE_RETVAL;
    let mut result = 0;

    for (i, element) in subgroup.iter().enumerate() {
        let found_log = value.equals(element);
        result |= (i as u8) & (found_log as u8);
        retval |= found_log;
    }

    (result, retval)
}
pub fn solve_dlp_order_power_of_five_explicit_subgroup<
    Fp2: Fp2Trait,
    const NUM_WORDS_ORDER: usize,
>(
    subgroup: &[Fp2; 5],
    value: &Fp2,
    e: usize,
    five_adic_basis: &[BigNum<NUM_WORDS_ORDER>],
) -> (BigNum<NUM_WORDS_ORDER>, u32) {
    let mut retval = SUCCESS_RETVAL;

    let prime_order_subgroup = subgroup.map(|element| {
        element.pow(
            &five_adic_basis[e - 1].to_le_bytes(),
            five_adic_basis[e - 1].nbits(),
        )
    });

    let mut partial_sum = BigNum::zero();
    for i in 0..e {
        let r = *value
            * subgroup[1]
                .pow(&partial_sum.to_le_bytes(), partial_sum.nbits())
                .invert(); // FIXME: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            &five_adic_basis[e - i - 1].to_le_bytes(),
            five_adic_basis[e - i - 1].nbits(),
        );

        let (x, ok) = solve_dlp_order_five_explicit_subgroup(&prime_order_subgroup, &u);
        partial_sum += (x as u64) * &five_adic_basis[i];
        retval &= ok;
    }

    (partial_sum, retval)
}

#[cfg(test)]
mod tests {
    use core::array;

    use rstest::rstest;

    use super::{SUCCESS_RETVAL, solve_dlp_small_prime_power_order};
    use crate::bn::BigNum;

    // p = 2 * 3^3 * 5^3 * 7^4 * 113 * 379 * 503 * 52837 + 1, such that |p| > 64 bits and Fp^2 has a multiplicative subgroup of order 5^3
    const P: [u64; 2] = [0x0000000000000c3f, 0x0000000000000001];
    fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = P,);

    #[rstest]
    fn test_simple_dlp() {
        let five_adic_basis: [BigNum<1>; 4] =
            array::from_fn(|exp| BigNum::from_prime_power(5, exp));

        // Generates the order-5^3 subgroup
        let generator = Fp2::const_decode_no_check(
            &[47, 47, 195, 76, 26, 17, 18, 144, 0],
            &[0; Fp::ENCODED_LENGTH],
        );
        // g^7
        let value = Fp2::const_decode_no_check(
            &[242, 200, 149, 124, 105, 197, 133, 202, 0],
            &[0; Fp::ENCODED_LENGTH],
        );

        assert_eq!(
            generator.pow(&[125], 7).equals(&Fp2::ONE),
            SUCCESS_RETVAL,
            "Generator doesn't have order 5^3"
        );

        let (x, ok) = solve_dlp_small_prime_power_order(&generator, &value, 5, 3, &five_adic_basis);

        assert_eq!(ok, SUCCESS_RETVAL);
        assert_eq!(x.to_le_bytes(), &[7]);
        assert_eq!(x.nbits(), 3);
    }

    #[rstest]
    fn test_complex_dlp() {
        let seven_adic_basis: [BigNum<1>; 5] =
            array::from_fn(|exp| BigNum::from_prime_power(7, exp));

        // Generates the order-7^4 subgroup
        let generator = Fp2::const_decode_no_check(
            &[82, 1, 230, 158, 56, 88, 5, 102, 0],
            &[0; Fp::ENCODED_LENGTH],
        );
        // g^501
        let value = Fp2::const_decode_no_check(
            &[148, 122, 116, 74, 152, 105, 50, 33, 0],
            &[0; Fp::ENCODED_LENGTH],
        );

        assert_eq!(
            generator.pow(&[97, 9], 12).equals(&Fp2::ONE),
            SUCCESS_RETVAL,
            "Generator doesn't have order 7^4"
        );

        let (x, ok) =
            solve_dlp_small_prime_power_order(&generator, &value, 7, 4, &seven_adic_basis);

        assert_eq!(ok, SUCCESS_RETVAL);
        assert_eq!(x.to_le_bytes(), &[255, 1]);
        assert_eq!(x.nbits(), 9);
    }
}
