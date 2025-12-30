use core::array;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{curve::Curve, projective_point::Point};
use num_bigint::{BigUint, RandBigInt as _};

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

pub fn sample_random_element_mod<const NUM_WORDS_MOD: usize>(
    modulus: &BigNum<NUM_WORDS_MOD>,
) -> BigNum<NUM_WORDS_MOD> {
    let mut rng = old_rand::thread_rng();

    let modulus = BigUint::from_bytes_le(&modulus.to_le_bytes());
    let element = rng.gen_biguint_below(&modulus);

    BigNum::new(&element.to_u64_digits())
    // BigNum::new(&[1])
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

    // let unit = BigNum::new(&[1]);
    BigNum::new(&unit.to_u64_digits())
}

// FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
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
    // TODO: is this valid? I would assume the operations between 3 random numbers also gives a random number. Prove this
    let cross_term = (&modulus - ((&matrix[0][1] * &matrix[1][0]) % &modulus)) % &modulus;
    let mut element = rng.gen_biguint_below(&modulus);
    while (&cross_term + (&matrix[0][0] * &element) % &modulus) % &modulus_base == BigUint::ZERO {
        element = rng.gen_biguint_below(&modulus);
    }
    matrix[1][1] = element;

    matrix.map(|row| row.map(|element| BigNum::new(&element.to_u64_digits())))
    // [
    //     [BigNum::new(&[1]), BigNum::zero()],
    //     [BigNum::zero(), BigNum::new(&[1])],
    // ]
}

// Randomly find a basis of the given torsion subgroup on the given curve
pub fn sample_random_torsion_basis_order_prime_power<
    Fp2: Fp2Trait,
    const NUM_WORDS_TORSION: usize,
    const NUM_WORDS_COF: usize,
>(
    curve: &Curve<Fp2>,
    torsion_subgroup_order_base: u8,
    torsion_subgroup_order: &BigNum<NUM_WORDS_TORSION>,
    order_cofactor: &BigNum<NUM_WORDS_COF>,
) -> (Point<Fp2>, Point<Fp2>, Fp2) {
    let mut rng = rand::rng();

    let torsion_subgroup_order_base = BigUint::from(torsion_subgroup_order_base);

    // TODO: include in paper WHY we can just divide p^e by p once to obtain a check that the point has exactly the order we need
    let reduced_torsion_subgroup_order =
        &BigUint::from_bytes_le(&torsion_subgroup_order.to_le_bytes())
            / &torsion_subgroup_order_base;
    let reduced_torsion_subgroup_order =
        BigNum::<NUM_WORDS_TORSION>::new(&reduced_torsion_subgroup_order.to_u64_digits());

    // Generate a point of the desired order
    // FIXME: should I break the loop condition only at the very end, all conditions being tested at once?
    // Or is this not even necessary, because exiting early only leaks information we know?
    // TODO: check how long it takes for these loops to terminate on average, as we are not generating deterministically as in the Sage implementation
    let U = loop {
        let U = curve.rand_point(&mut rng);

        /* Check point is in the torsion subgroup */

        let U_in_torsion_subgroup =
            curve.mul(&U, &order_cofactor.to_le_bytes(), order_cofactor.nbits());
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if U_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let U_saturated = curve.mul(
            &U_in_torsion_subgroup,
            &reduced_torsion_subgroup_order.to_le_bytes(),
            reduced_torsion_subgroup_order.nbits(),
        );
        if U_saturated.is_zero() == SUCCESS_RETVAL {
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

        let V_saturated = curve.mul(
            &V_in_torsion_subgroup,
            &reduced_torsion_subgroup_order.to_le_bytes(),
            reduced_torsion_subgroup_order.nbits(),
        );
        if V_saturated.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let V = V_in_torsion_subgroup;
        let UV = curve.sub(&U, &V);

        /* Check point is linearly independent to U */

        // FIXME: why the heck does this Weil pairing not produce the same as what Sage does??
        // In the hardcoded example below, it produces the square of what Sage does, for example
        let eUV = curve.weil_pairing(
            &U.to_pointx().x(),
            &V.to_pointx().x(),
            &UV.to_pointx().x(),
            &torsion_subgroup_order.to_le_bytes(),
            torsion_subgroup_order.nbits(),
        );

        let eUV_saturated = eUV.pow(
            &reduced_torsion_subgroup_order.to_le_bytes(),
            reduced_torsion_subgroup_order.nbits(),
        );
        if eUV_saturated.equals(&Fp2::ONE) == SUCCESS_RETVAL {
            continue;
        }
        // TODO: is this check necessary? Because of the fact that the group might have order m*n
        if eUV_saturated
            .pow(
                &torsion_subgroup_order_base.to_bytes_le(),
                torsion_subgroup_order_base.bits()
                    .try_into()
                    .expect("Size in bits of constant 5 is too big to fit in a usize (we do not ever expect this to happen)"),
            )
            .equals(&Fp2::ONE)
            == FAILURE_RETVAL
        {
            continue;
        }

        break (V, eUV);
    };

    (U, V, eUV)
}
