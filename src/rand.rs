use core::array;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{curve::Curve, projective_point::Point};
use num_bigint::{BigUint, RandBigInt as _};

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

pub fn sample_random_element_mod(modulus: &BigNum) -> BigNum {
    let mut rng = old_rand::thread_rng();

    let modulus = BigUint::from_bytes_le(modulus.as_le_bytes());
    let element = rng.gen_biguint_below(&modulus);

    BigNum::new(&element.to_u64_digits())
    // BigNum::new(&[1])
}

pub fn sample_random_unit_mod(modulus_base: u8, modulus: &BigNum) -> BigNum {
    let mut rng = old_rand::thread_rng();

    // Keep generating elements until we find an invertible one
    let modulus = BigUint::from_bytes_le(modulus.as_le_bytes());
    let mut unit = rng.gen_biguint_below(&modulus);
    while &unit % BigUint::from(modulus_base) == BigUint::ZERO {
        unit = rng.gen_biguint_below(&modulus);
    }

    let unit = BigNum::new(&unit.to_u64_digits());
    // let unit = BigNum::new(&[1]);

    unit
}

// FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
pub fn sample_random_invertible_matrix_mod(modulus_base: u8, modulus: &BigNum) -> [[BigNum; 2]; 2] {
    let mut rng = old_rand::thread_rng();

    let ONE = BigUint::from(1u8);
    let modulus_base = BigUint::from(modulus_base);
    let modulus = BigUint::from_bytes_le(modulus.as_le_bytes());

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
pub fn sample_random_torsion_basis<Fp2: Fp2Trait>(
    curve: &Curve<Fp2>,
    torsion_subgroup_order_base: u8,
    torsion_subgroup_order: &BigNum,
    order_cofactor: &BigNum,
) -> (Point<Fp2>, Point<Fp2>, Fp2) {
    let mut rng = rand::rng();

    let torsion_subgroup_order_base = BigUint::from(torsion_subgroup_order_base);

    // TODO: include in paper WHY we can just divide p^e by p once to obtain a check that the point has exactly the order we need
    let reduced_torsion_subgroup_order =
        &BigUint::from_bytes_le(torsion_subgroup_order.as_le_bytes())
            / &torsion_subgroup_order_base;
    let reduced_torsion_subgroup_order =
        BigNum::new(&reduced_torsion_subgroup_order.to_u64_digits());

    // Generate a point of the desired order
    // FIXME: should I break the loop condition only at the very end, all conditions being tested at once?
    // Or is this not even necessary, because exiting early only leaks information we know?
    // TODO: check how long it takes for these loops to terminate on average, as we are not generating deterministically as in the Sage implementation
    let U = loop {
        let U = curve.rand_point(&mut rng);

        /* Check point is in the torsion subgroup */

        let U_in_torsion_subgroup =
            curve.mul(&U, order_cofactor.as_le_bytes(), order_cofactor.nbits());
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if U_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let U_saturated = curve.mul(
            &U_in_torsion_subgroup,
            reduced_torsion_subgroup_order.as_le_bytes(),
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
            curve.mul(&V, order_cofactor.as_le_bytes(), order_cofactor.nbits());
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if V_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let V_saturated = curve.mul(
            &V_in_torsion_subgroup,
            reduced_torsion_subgroup_order.as_le_bytes(),
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
            torsion_subgroup_order.as_le_bytes(),
            torsion_subgroup_order.nbits(),
        );

        let eUV_saturated = eUV.pow(
            reduced_torsion_subgroup_order.as_le_bytes(),
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

    // let xU = PointX::from_x_coord(&Fp2::decode_reduce(&[
    //     3, 36, 23, 208, 115, 217, 172, 163, 98, 15, 91, 203, 150, 201, 77, 119, 81, 212, 9, 247,
    //     201, 65, 175, 184, 119, 101, 95, 190, 176, 162, 118, 222, 157, 140, 41, 11, 158, 65, 175,
    //     148, 126, 87, 155, 150, 173, 161, 130, 10, 199, 233, 219, 164, 247, 6, 60, 191, 212, 194,
    //     14, 110, 228, 248, 180, 81, 144, 208, 42, 227, 178, 92, 233, 102, 249, 149, 44, 248, 95,
    //     45, 181, 105, 58, 109, 45, 63, 187, 75, 102, 141, 143, 114, 62, 234, 107, 47, 204, 100,
    //     225, 156, 156, 112, 49, 219, 59, 73, 234, 241, 122, 2,
    // ]));
    // let xV = PointX::from_x_coord(&Fp2::decode_reduce(&[
    //     29, 148, 125, 128, 17, 251, 215, 243, 121, 68, 253, 137, 133, 48, 132, 193, 207, 95, 26,
    //     36, 39, 41, 163, 107, 95, 141, 213, 78, 17, 194, 51, 173, 239, 160, 184, 74, 81, 137, 37,
    //     148, 250, 220, 55, 41, 185, 4, 28, 108, 246, 213, 71, 72, 173, 55, 244, 213, 77, 72, 119,
    //     226, 154, 113, 30, 184, 78, 170, 0, 205, 112, 6, 254, 122, 217, 26, 224, 18, 135, 247, 104,
    //     59, 240, 58, 59, 233, 200, 84, 230, 241, 188, 95, 182, 145, 43, 206, 113, 116, 219, 162,
    //     177, 16, 238, 187, 54, 11, 143, 63, 170, 15,
    // ]));
    // let xUV = PointX::from_x_coord(&Fp2::decode_reduce(&[
    //     32, 44, 242, 20, 116, 222, 35, 180, 40, 3, 254, 145, 117, 21, 79, 35, 71, 51, 95, 179, 235,
    //     224, 244, 16, 106, 17, 43, 210, 41, 156, 219, 189, 91, 2, 155, 248, 10, 71, 75, 152, 43, 5,
    //     3, 83, 69, 232, 160, 140, 194, 199, 16, 199, 93, 103, 89, 35, 213, 13, 4, 62, 89, 157, 113,
    //     58, 44, 170, 0, 247, 163, 235, 153, 140, 182, 207, 136, 108, 160, 246, 139, 133, 99, 215,
    //     20, 210, 66, 23, 0, 66, 171, 119, 132, 162, 254, 163, 180, 213, 139, 177, 233, 249, 76,
    //     216, 73, 1, 64, 110, 254, 76,
    // ]));
    // let (U, V) = curve.lift_basis(&BasisX::from_points(&xU, &xV, &xUV));

    // // FIXME: why the heck does the Weil pairing produce the square of what Sage produces??
    // let eUV = curve.weil_pairing(
    //     &xU.x(),
    //     &xV.x(),
    //     &xUV.x(),
    //     torsion_subgroup_order.as_le_bytes(),
    //     torsion_subgroup_order.nbits(),
    // );
    // assert_eq!(
    //     eUV.pow(
    //         reduced_torsion_subgroup_order.as_le_bytes(),
    //         reduced_torsion_subgroup_order.nbits(),
    //     )
    //     .equals(&Fp2::ONE),
    //     FAILURE_RETVAL,
    //     "",
    // );
    // assert_eq!(
    //     eUV.pow(
    //         torsion_subgroup_order.as_le_bytes(),
    //         torsion_subgroup_order.nbits(),
    //     )
    //     .equals(&Fp2::ONE),
    //     SUCCESS_RETVAL,
    //     "",
    // );

    (U, V, eUV)
}
