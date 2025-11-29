use core::array;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{curve::Curve, projective_point::Point};
use num_bigint::{BigUint, RandBigInt as _};

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

pub fn sample_random_element_mod(modulus: &BigUint) -> BigNum {
    let mut rng = old_rand::thread_rng();

    let element = rng.gen_biguint_below(modulus);

    BigNum::new(&element.to_u64_digits())
}

pub fn sample_random_unit_mod(modulus: &BigUint) -> (BigNum, BigNum) {
    let mut rng = old_rand::thread_rng();

    // Keep generating elements until we find an invertible one
    let mut element = rng.gen_biguint_below(modulus);
    let mut inverse = element.modinv(modulus);
    while let None = inverse {
        element = rng.gen_biguint_below(modulus);
        inverse = element.modinv(modulus);
    }
    let Some(inverse) = inverse else {
        unreachable!("At this point, we are ensured to have an invertible element");
    };

    let element = BigNum::new(&element.to_u64_digits());

    let inverse = BigNum::new(&inverse.to_u64_digits());

    (element, inverse)
}

// FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
pub fn sample_random_invertible_matrix_mod(
    modulus_base: &BigUint,
    modulus: &BigUint,
) -> [[BigNum; 2]; 2] {
    let mut rng = old_rand::thread_rng();

    let ONE = BigUint::from(1u8);

    // Randomly generate the first 3 elements
    let mut matrix = array::from_fn(|_| [BigUint::ZERO; 2]);
    matrix[0][0] = rng.gen_biguint_range(&ONE, modulus); // This avoids getting a zero-term for the term in the determinant with our degree of freedom
    matrix[0][1] = rng.gen_biguint_below(modulus);
    matrix[1][0] = rng.gen_biguint_below(modulus);
    while ((&matrix[0][1] * &matrix[1][0]) % modulus_base) == BigUint::ZERO
        && ((modulus - &matrix[0][0]) % modulus_base) == BigUint::ZERO
    {
        matrix[0][0] = rng.gen_biguint_range(&ONE, modulus); // This avoids getting a zero-term for the term in the determinant with our degree of freedom
    }

    // Select the 4th element to have gcd(det(D), 5^c) == 1
    // TODO: is this valid? I would assume the operations between 3 random numbers also gives a random number. Prove this
    let cross_term = (modulus - ((&matrix[0][1] * &matrix[1][0]) % modulus)) % modulus;
    let mut element = rng.gen_biguint_below(modulus);
    let mut det_inverse =
        ((&cross_term + (&matrix[0][0] * &element) % modulus) % modulus).modinv(modulus);
    while let None = det_inverse {
        element = rng.gen_biguint_below(modulus);
        det_inverse =
            ((&cross_term + (&matrix[0][0] * &element) % modulus) % modulus).modinv(modulus);
    }
    matrix[1][1] = element;

    matrix.map(|row| row.map(|element| BigNum::new(&element.to_u64_digits())))
}

// Randomly find a basis of the given torsion subgroup on the given curve
pub fn sample_random_torsion_basis<Fp2: Fp2Trait>(
    curve: &Curve<Fp2>,
    torsion_subgroup_order: &BigUint,
    order_cofactor: &BigUint,
) -> (Point<Fp2>, Point<Fp2>, Fp2) {
    let mut rng = rand::rng();

    let FIVE = BigUint::from(5u8);

    let reduced_torsion_subgroup_order = torsion_subgroup_order / &FIVE;
    let reduced_torsion_subgroup_order_bitsize = reduced_torsion_subgroup_order
        .bits()
        .try_into()
        .expect("Size in bits of the (torsion subgroup order / 5) is too big to fit in a usize (we do not ever expect this to happen)");
    let reduced_torsion_subgroup_order = reduced_torsion_subgroup_order.to_bytes_le();
    let torsion_subgroup_order_bitsize = torsion_subgroup_order
        .bits()
        .try_into()
        .expect("Size in bits of the torsion subgroup order is too big to fit in a usize (we do not ever expect this to happen)");
    let torsion_subgroup_order = torsion_subgroup_order.to_bytes_le();
    let order_cofactor_bitsize = order_cofactor
        .bits()
        .try_into()
        .expect("Size in bits of the cofactor (p + 1)/(pi)^(ei) is too big to fit in a usize (we do not ever expect this to happen)");
    let order_cofactor = order_cofactor.to_bytes_le();

    // Generate a point of the desired order
    // FIXME: should I break the loop condition only at the very end, all conditions being tested at once?
    // Or is this not even necessary, because exiting early only leaks information we know?
    // TODO: check how long it takes for these loops to terminate on average, as we are not generating deterministically as in the Sage implementation
    let U = loop {
        let U = curve.rand_point(&mut rng);

        /* Check point is in the torsion subgroup */

        let U_in_torsion_subgroup = curve.mul(&U, &order_cofactor, order_cofactor_bitsize);
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if U_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let U_saturated = curve.mul(
            &U_in_torsion_subgroup,
            &reduced_torsion_subgroup_order,
            reduced_torsion_subgroup_order_bitsize,
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

        let V_in_torsion_subgroup = curve.mul(&V, &order_cofactor, order_cofactor_bitsize);
        // We don't want a point in the ((p + 1)/(pi)^(ei))-torsion
        if V_in_torsion_subgroup.is_zero() == SUCCESS_RETVAL {
            continue;
        }

        let V_saturated = curve.mul(
            &V_in_torsion_subgroup,
            &reduced_torsion_subgroup_order,
            reduced_torsion_subgroup_order_bitsize,
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
            &torsion_subgroup_order,
            torsion_subgroup_order_bitsize,
        );

        let eUV_saturated = eUV.pow(
            &reduced_torsion_subgroup_order,
            reduced_torsion_subgroup_order_bitsize,
        );
        if eUV_saturated.equals(&Fp2::ONE) == SUCCESS_RETVAL {
            continue;
        }
        // TODO: is this check necessary? Because of the fact that the group might have order m*n
        if eUV_saturated
            .pow(
                &FIVE.to_bytes_le(),
                FIVE.bits()
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
    //     74, 197, 179, 60, 61, 71, 181, 238, 134, 215, 70, 233, 42, 125, 178, 160, 122, 74, 93, 70,
    //     55, 174, 167, 240, 208, 192, 66, 113, 92, 61, 229, 3, 17, 243, 39, 220, 41, 221, 102, 217,
    //     166, 69, 165, 197, 137, 6, 178, 108, 68, 208, 75, 132, 29, 98, 211, 72, 56, 202, 64, 195,
    //     8, 61, 46, 181, 95, 190, 226, 231, 4, 46, 74, 83, 140, 103, 164, 178, 47, 92, 94, 189, 71,
    //     117, 202, 53, 190, 254, 253, 179, 113, 253, 146, 68, 86, 121, 139, 209, 76, 219, 209, 231,
    //     93, 207, 56, 147, 13, 176, 187, 23,
    // ]));
    // let xV = PointX::from_x_coord(&Fp2::decode_reduce(&[
    //     135, 58, 149, 162, 71, 21, 110, 65, 140, 176, 190, 97, 35, 153, 221, 164, 27, 174, 240,
    //     122, 245, 113, 63, 36, 40, 23, 23, 140, 171, 173, 142, 244, 239, 208, 180, 220, 186, 145,
    //     41, 170, 15, 153, 170, 69, 252, 133, 7, 40, 6, 45, 87, 57, 210, 53, 64, 102, 36, 72, 20,
    //     46, 107, 111, 45, 167, 80, 189, 188, 94, 3, 86, 165, 82, 65, 248, 166, 14, 216, 110, 83,
    //     136, 51, 229, 88, 169, 189, 212, 120, 246, 217, 164, 43, 110, 126, 94, 137, 136, 190, 101,
    //     147, 171, 33, 209, 47, 152, 120, 48, 24, 4,
    // ]));
    // let xUV = PointX::from_x_coord(&Fp2::decode_reduce(&[
    //     77, 215, 149, 22, 64, 114, 182, 118, 22, 44, 208, 178, 202, 144, 130, 204, 125, 146, 126,
    //     148, 179, 17, 148, 187, 194, 164, 234, 12, 211, 89, 111, 179, 29, 140, 77, 29, 122, 199,
    //     151, 202, 9, 61, 167, 8, 163, 96, 2, 181, 233, 40, 155, 155, 140, 11, 119, 224, 139, 100,
    //     109, 18, 235, 101, 77, 160, 132, 141, 177, 236, 26, 195, 85, 101, 232, 74, 106, 250, 131,
    //     167, 196, 87, 230, 73, 213, 255, 88, 44, 47, 45, 124, 185, 134, 97, 87, 55, 160, 114, 182,
    //     166, 9, 49, 202, 232, 100, 244, 139, 212, 204, 57,
    // ]));
    // let (U, V) = curve.lift_basis(&BasisX::from_points(&xU, &xV, &xUV));

    // // FIXME: why the heck does the Weil pairing produce the square of what Sage produces??
    // let (eUV, ok) = curve
    //     .weil_pairing(
    //         &xU.x(),
    //         &xV.x(),
    //         &xUV.x(),
    //         &torsion_subgroup_order,
    //         torsion_subgroup_order_bitsize,
    //     )
    //     .sqrt();
    // assert_eq!(ok, SUCCESS_RETVAL, "Weil pairing doesn't produce a square");
    // assert_eq!(
    //     eUV.pow(
    //         &reduced_torsion_subgroup_order,
    //         reduced_torsion_subgroup_order_bitsize
    //     )
    //     .equals(&Fp2::ONE),
    //     FAILURE_RETVAL
    // );
    // assert_eq!(
    //     eUV.pow(&torsion_subgroup_order, torsion_subgroup_order_bitsize)
    //         .equals(&Fp2::ONE),
    //     SUCCESS_RETVAL
    // );

    (U, V, eUV)
}
