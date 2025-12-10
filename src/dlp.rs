use fp2::traits::Fp2 as Fp2Trait;

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

pub fn solve_dlp_small_prime_order<Fp2: Fp2Trait>(
    generator: &Fp2,
    value: &Fp2,
    p: usize,
) -> (usize, u32) {
    let mut retval = FAILURE_RETVAL;

    let mut element = Fp2::ONE;
    let mut result = 0;
    for i in 0..p {
        let found_log = value.equals(&element);
        result |= i & (((found_log as usize) << 32) | found_log as usize);
        retval |= found_log;

        element *= *generator;
    }

    (result, retval)
}

// WARN: Vartime with respect to public parameters (the order of the 5^c-torsion subgroup)
pub fn solve_dlp_small_prime_power_order<Fp2: Fp2Trait>(
    generator: &Fp2,
    value: &Fp2,
    p: usize,
    e: usize,
) -> (BigNum, u32) {
    let mut retval = SUCCESS_RETVAL;

    let p_to_the_e_basis = (0..=e)
        .map(|exp| BigNum::from_prime_power(p, exp))
        .collect::<Vec<_>>();

    let prime_order_subgroup_generator = generator.pow(
        p_to_the_e_basis[e - 1].as_le_bytes(),
        p_to_the_e_basis[e - 1].nbits(),
    );
    assert_eq!(
        prime_order_subgroup_generator.equals(&Fp2::ONE),
        FAILURE_RETVAL,
        "g has order < {}^{}",
        p,
        e,
    );

    let mut partial_solutions = Vec::with_capacity(e);
    let mut partial_sum = BigNum::zero();
    for i in 0..e {
        let r = *value
            * generator
                .pow(partial_sum.as_le_bytes(), partial_sum.nbits())
                .invert(); // TODO: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            p_to_the_e_basis[e - i - 1].as_le_bytes(),
            p_to_the_e_basis[e - i - 1].nbits(),
        );

        let (x, ok) = solve_dlp_small_prime_order(&prime_order_subgroup_generator, &u, p);
        partial_solutions.push(x);
        partial_sum = &partial_sum + &((x as u64) * &p_to_the_e_basis[i]);
        retval &= ok;
    }

    (partial_sum, retval)
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

    element2.set_invert();

    let found_log = value.equals(&element2);
    result |= 3 & (found_log as u8);
    retval |= found_log;

    element.set_invert();

    let found_log = value.equals(&element);
    result |= 4 & (found_log as u8);
    retval |= found_log;

    (result, retval)
}

// WARN: Vartime with respect to public parameters (the order of the 5^c-torsion subgroup)
pub fn solve_dlp_order_power_of_five<Fp2: Fp2Trait>(
    generator: &Fp2,
    value: &Fp2,
    e: usize,
) -> (BigNum, u32) {
    let mut retval = SUCCESS_RETVAL;

    let p_to_the_e_basis = (0..=e)
        .map(|exp| BigNum::from_prime_power(5, exp))
        .collect::<Vec<_>>();

    let prime_order_subgroup_generator = generator.pow(
        p_to_the_e_basis[e - 1].as_le_bytes(),
        p_to_the_e_basis[e - 1].nbits(),
    );
    assert_eq!(
        prime_order_subgroup_generator.equals(&Fp2::ONE),
        FAILURE_RETVAL,
        "g has order < 5^{}",
        e,
    );

    let mut partial_solutions = Vec::with_capacity(e);
    let mut partial_sum = BigNum::zero();
    for i in 0..e {
        let r = *value
            * generator
                .pow(partial_sum.as_le_bytes(), partial_sum.nbits())
                .invert(); // TODO: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            p_to_the_e_basis[e - i - 1].as_le_bytes(),
            p_to_the_e_basis[e - i - 1].nbits(),
        );

        let (x, ok) = solve_dlp_order_five(&prime_order_subgroup_generator, &u);
        partial_solutions.push(x);
        partial_sum = &partial_sum + &((x as u64) * &p_to_the_e_basis[i]);
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
pub fn solve_dlp_order_power_of_five_explicit_subgroup<Fp2: Fp2Trait>(
    subgroup: &[Fp2; 5],
    value: &Fp2,
    e: usize,
) -> (BigNum, u32) {
    let mut retval = SUCCESS_RETVAL;

    let p_to_the_e_basis = (0..=e)
        .map(|exp| BigNum::from_prime_power(5, exp))
        .collect::<Vec<_>>();

    let prime_order_subgroup = subgroup.map(|element| {
        element.pow(
            p_to_the_e_basis[e - 1].as_le_bytes(),
            p_to_the_e_basis[e - 1].nbits(),
        )
    });

    let mut partial_solutions = Vec::with_capacity(e);
    let mut partial_sum = BigNum::zero();
    for i in 0..e {
        let r = *value
            * subgroup[1]
                .pow(partial_sum.as_le_bytes(), partial_sum.nbits())
                .invert(); // TODO: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            p_to_the_e_basis[e - i - 1].as_le_bytes(),
            p_to_the_e_basis[e - i - 1].nbits(),
        );

        partial_solutions.push(x);
        let (x, ok) = solve_dlp_order_five_explicit_subgroup(&prime_order_subgroup, &u);
        partial_sum = &partial_sum + &((x as u64) * &p_to_the_e_basis[i]);
        retval &= ok;
    }

    (partial_sum, retval)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::{SUCCESS_RETVAL, solve_dlp_small_prime_power_order};

    // p = 2 * 3^3 * 5^3 * 7^4 * 113 * 379 * 503 * 52837 + 1, such that |p| > 64 bits and Fp^2 has a multiplicative subgroup of order 5^3
    const P: [u64; 2] = [0x0000000000000c3f, 0x0000000000000001];
    fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = P,);

    #[rstest]
    fn test_simple_dlp() {
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

        let (x, ok) = solve_dlp_small_prime_power_order(&generator, &value, 5, 3);

        assert_eq!(ok, SUCCESS_RETVAL);
        assert_eq!(x.as_le_bytes(), &[7]);
        assert_eq!(x.nbits(), 3);
    }
    #[rstest]
    fn test_complex_dlp() {
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

        let (x, ok) = solve_dlp_small_prime_power_order(&generator, &value, 7, 4);

        assert_eq!(ok, SUCCESS_RETVAL);
        assert_eq!(x.as_le_bytes(), &[255, 1]);
        assert_eq!(x.nbits(), 9);
    }
}
