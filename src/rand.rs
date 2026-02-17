use core::array;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{curve::Curve, projective_point::Point};
use rand::Rng as _;

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

pub fn sample_random_element_mod<const NUM_WORDS_MOD: usize>(
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> BigNum<NUM_WORDS_MOD> {
    let mut rng = rand::rng();

    let actual_num_words = modulus.nbits().div_ceil(u64::BITS as usize);
    let zero_out_mask =
        u64::MAX >> (u64::BITS - ((modulus.nbits() - 1) % (u64::BITS as usize)) as u32 - 1);

    let mut words = [0; NUM_WORDS_MOD];
    rng.fill(&mut words[..actual_num_words]);
    words[actual_num_words - 1] &= zero_out_mask;

    let mut element = BigNum::new(&words);
    while &element >= modulus {
        rng.fill(&mut words[..actual_num_words]);
        words[actual_num_words - 1] &= zero_out_mask;

        element = BigNum::new(&words);
    }

    element
}

pub fn sample_random_unit_mod_power_of_two<const NUM_WORDS_MOD: usize>(
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> BigNum<NUM_WORDS_MOD> {
    let mut unit = sample_random_element_mod(modulus);
    while unit.is_divisible_by_two() {
        unit = sample_random_element_mod(modulus);
    }

    unit
}

// WARN: modulus_base must satisfy 2^64 == 1 (mod modulus_base)
pub fn sample_random_unit_mod_special_prime_power<const NUM_WORDS_MOD: usize>(
    modulus_base: u8,
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> BigNum<NUM_WORDS_MOD> {
    let mut unit = sample_random_element_mod(modulus);
    while unit.is_divisible_by_special_prime(modulus_base) {
        unit = sample_random_element_mod(modulus);
    }

    unit
}

// WARN: All primes in special_primes_coprime_to must satisfy 2^64 == 1 (mod p)
pub fn sample_random_secret_degree<const NUM_WORDS_2: usize>(
    effective_two_torsion: &BigNum<NUM_WORDS_2>,
    special_primes_coprime_to: &[u8],
) -> (BigNum<NUM_WORDS_2>, BigNum<NUM_WORDS_2>) {
    // Keep generating elements until we find one such that q*(2^a - q) is coprime to all given primes
    let mut element = sample_random_unit_mod_power_of_two(effective_two_torsion);
    let mut dual_element = effective_two_torsion - &element;
    while special_primes_coprime_to.iter().any(|&p| {
        element.is_divisible_by_special_prime(p) || dual_element.is_divisible_by_special_prime(p)
    }) {
        element = sample_random_unit_mod_power_of_two(&effective_two_torsion);
        dual_element = effective_two_torsion - &element;
    }

    (element, dual_element)
}

// FIXME: Validate that this algorithm generates a uniformly random invertible matrices in SL_2(Z_(5^c))
// WARN: modulus_base must satisfy 2^64 == 1 (mod modulus_base)
pub fn sample_random_invertible_matrix_mod_special_prime_power<const NUM_WORDS_MOD: usize>(
    modulus_base: u8,
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> [[BigNum<NUM_WORDS_MOD>; 2]; 2] {
    // Randomly generate the first 3 elements
    let mut matrix: [[BigNum<NUM_WORDS_MOD>; 2]; 2] =
        array::from_fn(|_| array::from_fn(|_| BigNum::<NUM_WORDS_MOD>::zero()));

    matrix[0][1] = sample_random_element_mod(modulus);
    matrix[1][0] = sample_random_element_mod(modulus);
    let cross_term = modulus - matrix[0][1].mul_mod(&matrix[1][0], modulus);

    // Select the 1st and 4th elements so that gcd(det(D), modulus) == 1, i.e. so D is invertible in modulus
    // FIXME: Is this the fastest way to do this? If we fix 3/4 of the elements with the necessary properties
    // to ensure success instead, and only regenerate the 4th, how many iterations do we save? What is the
    // expected number of iterations in both cases?
    matrix[0][0] = sample_random_element_mod(modulus);
    matrix[1][1] = sample_random_element_mod(modulus);
    while (matrix[0][0].mul_mod(&matrix[1][1], modulus) + &cross_term)
        .is_divisible_by_special_prime(modulus_base)
    {
        matrix[0][0] = sample_random_element_mod(modulus);
        matrix[1][1] = sample_random_element_mod(modulus);
    }

    matrix
}

// Randomly find a basis of the given torsion subgroup on the given curve
// WARN: Only works for torsion subgroup orders with prime factors < 256
pub fn sample_random_torsion_basis<
    Fp2: Fp2Trait,
    const NUM_PRIME_FACTORS: usize,
    const NUM_WORDS_TORSION: usize,
    const NUM_WORDS_COF: usize,
>(
    curve: &Curve<Fp2>,
    torsion_subgroup_order_prime_factors: &[u8; NUM_PRIME_FACTORS],
    torsion_subgroup_reduced_orders: &[&BigNum<NUM_WORDS_TORSION>; NUM_PRIME_FACTORS],
    torsion_subgroup_order: &BigNum<NUM_WORDS_TORSION>,
    order_cofactor: &BigNum<NUM_WORDS_COF>,
) -> (Point<Fp2>, Point<Fp2>, Fp2) {
    let mut rng = rand::rng();

    let torsion_subgroup_order_prime_factors =
        torsion_subgroup_order_prime_factors.map(|prime| BigNum::<1>::from_prime(prime.into()));


    // Generate a point of the desired order
    // ASK: Should I break the loop condition only at the very end, all conditions being tested at once?
    // Or is this not even necessary, because exiting early only leaks information we know?
    // My guess is that when doing rejection sampling, it doesn't matter when we reject, since we know
    // the final result must not satisfy any of the rejection criteria anyway.
    // TODO: Check how long it takes for these loops to terminate on average, as we are not generating deterministically as in the Sage implementation
    let U = loop {
        let U = curve.rand_point(&mut rng);

        /* Check point is in the torsion subgroup */

        let U_in_torsion_subgroup =
            curve.mul(&U, &order_cofactor.to_le_bytes(), order_cofactor.nbits());
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if U_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let mut Us_saturated = torsion_subgroup_reduced_orders.iter().map(|reduced_order| {
            curve.mul(
                &U_in_torsion_subgroup,
                &reduced_order.to_le_bytes(),
                reduced_order.nbits(),
            )
        });
        // ASK: Is it a problem that any() short-circuits? See above comment
        if Us_saturated.any(|U_saturated| U_saturated.is_zero() == SUCCESS_RETVAL) {
            continue;
        }

        break U_in_torsion_subgroup;
    };

    // Generate another point of the desired order, linearly independent to the first
    // Includes Weil pairing between both points
    let (V, eUV) = loop {
        let V = curve.rand_point(&mut rng);

        /* Check point is in the torsion subgroup */

        let V_in_torsion_subgroup =
            curve.mul(&V, &order_cofactor.to_le_bytes(), order_cofactor.nbits());
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if V_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let mut Vs_saturated = torsion_subgroup_reduced_orders.iter().map(|reduced_order| {
            curve.mul(
                &V_in_torsion_subgroup,
                &reduced_order.to_le_bytes(),
                reduced_order.nbits(),
            )
        });
        // ASK: Is it a problem that any() short-circuits? See above comment
        if Vs_saturated.any(|V_saturated| V_saturated.is_zero() == SUCCESS_RETVAL) {
            continue;
        }

        let V = V_in_torsion_subgroup;
        let UV = curve.sub(&U, &V);

        /* Check point is linearly independent to U */

        let eUV = curve.weil_pairing(
            &U.to_pointx().x(),
            &V.to_pointx().x(),
            &UV.to_pointx().x(),
            &torsion_subgroup_order.to_le_bytes(),
            torsion_subgroup_order.nbits(),
        );

        // ASK: Is it a problem that any() short-circuits? See above comment
        let eUVs_saturated = torsion_subgroup_reduced_orders
            .iter()
            .map(|reduced_order| eUV.pow(&reduced_order.to_le_bytes(), reduced_order.nbits()))
            .collect::<Vec<_>>();
        if eUVs_saturated
            .iter()
            .any(|eUV_saturated| eUV_saturated.equals(&Fp2::ONE) == SUCCESS_RETVAL)
        {
            continue;
        }
        // FIXME: Is this check necessary? Because of the fact that the group might have order m*n
        // ASK: Is it a problem that any() short-circuits? See above comment
        if eUVs_saturated
            .iter()
            .zip(torsion_subgroup_order_prime_factors.iter())
            .any(|(eUV_saturated, prime)| {
                eUV_saturated
                    .pow(&prime.to_le_bytes(), prime.nbits())
                    .equals(&Fp2::ONE)
                    == FAILURE_RETVAL
            })
        {
            continue;
        }

        break (V, eUV);
    };

    (U, V, eUV)
}
