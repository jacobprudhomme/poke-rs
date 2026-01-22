use core::array;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{curve::Curve, projective_point::Point};
use num_bigint::{BigUint, RandBigInt as _};

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

pub fn sample_random_secret_degree<const NUM_WORDS_2: usize>(
    effective_two_torsion: &BigNum<NUM_WORDS_2>,
    odd_primes_coprime_to: &[u8],
) -> (BigNum<NUM_WORDS_2>, BigNum<NUM_WORDS_2>) {
    let mut rng = old_rand::thread_rng();

    let TWO = BigUint::from(2u8);
    let coprime_to = odd_primes_coprime_to
        .iter()
        .map(|&n| BigUint::from(n))
        .collect::<Vec<_>>();

    // Keep generating elements until we find one such that q*(2^a - q) is coprime to all given primes
    let effective_two_torsion = BigUint::from_bytes_le(&effective_two_torsion.to_le_bytes());
    let mut element = rng.gen_biguint_below(&effective_two_torsion);
    let mut dual_element = &effective_two_torsion - &element;
    while &element % &TWO == BigUint::ZERO
        || coprime_to
            .iter()
            .any(|p| &element % p == BigUint::ZERO || &dual_element % p == BigUint::ZERO)
    {
        element = rng.gen_biguint_below(&effective_two_torsion);
        dual_element = &effective_two_torsion - &element;
    }

    (
        BigNum::new(&element.to_u64_digits()),
        BigNum::new(&dual_element.to_u64_digits()),
    )
}

pub fn sample_random_element_mod<const NUM_WORDS_MOD: usize>(
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> BigNum<NUM_WORDS_MOD> {
    let mut rng = old_rand::thread_rng();

    let modulus = BigUint::from_bytes_le(&modulus.to_le_bytes());
    let element = rng.gen_biguint_below(&modulus);

    BigNum::new(&element.to_u64_digits())
}

pub fn sample_random_unit_mod_prime_power<const NUM_WORDS_MOD: usize>(
    modulus_base: u8,
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> BigNum<NUM_WORDS_MOD> {
    let mut rng = old_rand::thread_rng();

    // Keep generating elements until we find an invertible one
    let modulus = BigUint::from_bytes_le(&modulus.to_le_bytes());
    let mut unit = rng.gen_biguint_below(&modulus);
    while &unit % BigUint::from(modulus_base) == BigUint::ZERO {
        unit = rng.gen_biguint_below(&modulus);
    }

    BigNum::new(&unit.to_u64_digits())
}

// FIXME: Validate that this algorithm generates a uniformly random invertible matrices in SL_2(Z_(5^c))
pub fn sample_random_invertible_matrix_mod_prime_power<const NUM_WORDS_MOD: usize>(
    modulus_base: u8,
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> [[BigNum<NUM_WORDS_MOD>; 2]; 2] {
    let mut rng = old_rand::thread_rng();

    let ONE = BigUint::from(1u8);
    let modulus_base = BigUint::from(modulus_base);
    let modulus = BigUint::from_bytes_le(&modulus.to_le_bytes());

    // Randomly generate the first 3 elements
    let mut matrix: [[BigUint; 2]; 2] = array::from_fn(|_| [BigUint::ZERO; 2]);
    matrix[0][0] = rng.gen_biguint_range(&ONE, &modulus); // This avoids getting a zero-term for the term in the determinant with our degree of freedom
    matrix[0][1] = rng.gen_biguint_below(&modulus);
    matrix[1][0] = rng.gen_biguint_below(&modulus);
    while ((&matrix[0][1] * &matrix[1][0]) % &modulus_base) == BigUint::ZERO
        && ((&modulus - &matrix[0][0]) % &modulus_base) == BigUint::ZERO
    {
        matrix[0][0] = rng.gen_biguint_range(&ONE, &modulus); // This avoids getting a zero-term for the term in the determinant with our degree of freedom
    }

    // Select the 4th element to have gcd(det(D), 5^c) == 1
    // FIXME: Is this valid? I would assume the operations between 3 random numbers also gives a random number. Prove this
    let cross_term = (&modulus - ((&matrix[0][1] * &matrix[1][0]) % &modulus)) % &modulus;
    let mut element = rng.gen_biguint_below(&modulus);
    while (&cross_term + (&matrix[0][0] * &element) % &modulus) % &modulus_base == BigUint::ZERO {
        element = rng.gen_biguint_below(&modulus);
    }
    matrix[1][1] = element;

    matrix.map(|row| row.map(|element| BigNum::new(&element.to_u64_digits())))
}

// Randomly find a basis of the given torsion subgroup on the given curve
// WARN: Only works for torsion subgroup orders with prime factors < 256
pub fn sample_random_torsion_basis<
    Fp2: Fp2Trait,
    const NUM_WORDS_TORSION: usize,
    const NUM_WORDS_COF: usize,
>(
    curve: &Curve<Fp2>,
    torsion_subgroup_order_prime_factors: &[u8],
    torsion_subgroup_order: &BigNum<NUM_WORDS_TORSION>,
    order_cofactor: &BigNum<NUM_WORDS_COF>,
) -> (Point<Fp2>, Point<Fp2>, Fp2) {
    let mut rng = rand::rng();

    let torsion_subgroup_order_prime_factors = torsion_subgroup_order_prime_factors
        .iter()
        .map(|&prime| BigUint::from(prime))
        .collect::<Vec<_>>();
    let torsion_subgroup_order_biguint =
        BigUint::from_bytes_le(&torsion_subgroup_order.to_le_bytes());

    let reduced_torsion_subgroup_orders = &torsion_subgroup_order_prime_factors
        .iter()
        .map(|prime| &torsion_subgroup_order_biguint / prime)
        .map(|reduced_order| BigNum::<NUM_WORDS_TORSION>::new(&reduced_order.to_u64_digits()))
        .collect::<Vec<_>>();

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

        let mut Us_saturated = reduced_torsion_subgroup_orders.iter().map(|reduced_order| {
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

        let mut Vs_saturated = reduced_torsion_subgroup_orders.iter().map(|reduced_order| {
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
        let eUVs_saturated = reduced_torsion_subgroup_orders
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
        if eUVs_saturated.iter().zip(torsion_subgroup_order_prime_factors.iter()).any(|(eUV_saturated, prime)| {
            eUV_saturated
                .pow(
                    &prime.to_bytes_le(),
                    prime.bits()
                        .try_into()
                        .expect("Size in bits of constant 5 is too big to fit in a usize (we do not ever expect this to happen)"),
                )
                .equals(&Fp2::ONE)
                == FAILURE_RETVAL
        }) {
            continue;
        }

        break (V, eUV);
    };

    (U, V, eUV)
}
