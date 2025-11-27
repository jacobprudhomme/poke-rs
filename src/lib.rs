#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]

use std::{array, io::Read as _, marker::PhantomData};

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
    utilities::bn::{
        add_bn_vartime, bn_bit_length_vartime, mul_bn_by_u64_vartime, prime_power_to_bn_vartime,
    },
};
use num_bigint::{BigUint, RandBigInt as _};
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _},
};

pub mod example_keypairs;
pub mod fields;
pub mod params;

pub const SUCCESS_RETVAL: u32 = u32::MAX;
pub const FAILURE_RETVAL: u32 = u32::MIN;

#[derive(Debug, PartialEq)]
pub struct BigNum {
    pub repr: Vec<u8>,
    pub bitlen: usize,
}

pub struct PublicParams<Fp2: Fp2Trait> {
    pub starting_curve: Curve<Fp2>,
    pub full_two_torsion_order: BigUint,
    pub full_two_torsion_exp: usize,
    pub effective_two_torsion_order: BigUint,
    pub effective_two_torsion_exp: usize,
    pub three_torsion_order: BigUint,
    pub three_torsion_exp: usize,
    pub five_torsion_order: BigUint,
    pub five_torsion_exp: usize,
    pub five_torsion_cofactor: BigUint,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
    pub five_torsion_basis: BasisX<Fp2>,
}

pub struct InkePublicParams<Fp2: Fp2Trait> {
    pub starting_curve: Curve<Fp2>,
    pub effective_two_torsion_order: BigUint,
    pub effective_two_torsion_exp: usize,
    pub three_torsion_order: BigUint,
    pub three_torsion_exp: usize,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
}

// FIXME: represent scalars as their LE byte arrays and bitsize. Removes external dependency on num-bigint
pub struct PrvKey<Fp2: Fp2Trait> {
    pub q: BigUint,
    pub alpha: BigUint,
    pub beta: BigUint,
    pub delta: BigUint,
    pub _field: PhantomData<Fp2>,
}

pub struct InkePrvKey<Fp2: Fp2Trait> {
    pub q: BigUint,
    pub alpha: BigUint,
    pub beta: BigUint,
    pub _field: PhantomData<Fp2>,
}

pub struct PubKey<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_img: BasisX<Fp2>,
    pub masked_three_torsion_basis_img: BasisX<Fp2>,
    pub masked_five_torsion_basis_img: BasisX<Fp2>,
}

pub struct InkePubKey<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_img: BasisX<Fp2>,
    pub masked_three_torsion_basis_img: BasisX<Fp2>,
    pub intermediate_curve: Curve<Fp2>,
    pub masked_three_torsion_basis_img_intermediate: BasisX<Fp2>,
}

pub struct Ciphertext<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EB: BasisX<Fp2>,
    pub masked_five_torsion_basis_EB: BasisX<Fp2>,
    pub shared_end_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EAB: BasisX<Fp2>,
    pub encrypted_message: Vec<u8>,
}

pub struct InkeCiphertext<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EB: BasisX<Fp2>,
    pub shared_end_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EAB: BasisX<Fp2>,
    pub encrypted_message: Vec<u8>,
}

pub fn byte_bn_to_word_bn(bytes: &BigNum) -> Vec<u64> {
    bytes
        .repr
        .chunks(8)
        .map(|word_bytes| {
            let mut word_bytes = word_bytes.to_vec();
            word_bytes.resize(8, 0);
            u64::from_le_bytes(word_bytes.try_into().unwrap())
        })
        .collect()
}

pub fn word_bn_to_byte_bn(words: &[u64]) -> BigNum {
    let bitlen = bn_bit_length_vartime(words);
    let mut bytes = words
        .iter()
        .flat_map(|word| word.to_le_bytes())
        .collect::<Vec<_>>();
    while bytes.len() > 1
        && let Some(&last_byte) = bytes.last()
        && last_byte == 0
    {
        bytes.pop();
    }

    BigNum {
        repr: bytes,
        bitlen,
    }
}

// WARN: Assumes we are given an a priori invertible element as input
fn invert_element_mod(element: &BigUint, modulus: &BigUint) -> BigNum {
    let inverse = element.modinv(modulus);
    let Some(inverse) = inverse else {
        unreachable!("We expect an invertible element as input");
    };

    let bitlen = inverse
        .bits()
        .try_into()
        .expect("Size in bits of the inverse scalar is too big to fit in a usize (we do not ever expect this to happen)");
    let repr = inverse.to_bytes_le();

    BigNum { repr, bitlen }
}

fn sample_random_element_mod(modulus: &BigUint) -> BigNum {
    let mut rng = old_rand::thread_rng();

    let element = rng.gen_biguint_below(modulus);

    // Transform element to our own BigNum type
    let bitlen = element
        .bits()
        .try_into()
        .expect("Size in bits of the scalar is too big to fit in a usize (we do not ever expect this to happen)");
    let repr = element.to_bytes_le();

    BigNum { repr, bitlen }
}

fn sample_random_unit_mod(modulus: &BigUint) -> (BigNum, BigNum) {
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

    // Transform element to our own BigNum type
    let element_bitlen = element
        .bits()
        .try_into()
        .expect("Size in bits of the scalar is too big to fit in a usize (we do not ever expect this to happen)");
    let element = element.to_bytes_le();
    let element = BigNum {
        repr: element,
        bitlen: element_bitlen,
    };

    // Transform inverse of element to our own BigNum type
    let inverse_bitlen = inverse
        .bits()
        .try_into()
        .expect("Size in bits of the inverse scalar is too big to fit in a usize (we do not ever expect this to happen)");
    let inverse = inverse.to_bytes_le();
    let inverse = BigNum {
        repr: inverse,
        bitlen: inverse_bitlen,
    };

    (element, inverse)
}

// FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
fn sample_random_invertible_matrix_mod(
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

    matrix.map(|row| {
        row.map(|element| {
            let bitlen = element
                .bits()
                .try_into()
                .expect("Size in bits of the matrix element is too big to fit in a usize (we do not ever expect this to happen)");
            let repr = element.to_bytes_le();

            BigNum { repr, bitlen }
        })
    })
}

pub fn solve_dlp_small_prime_order<Fp2: Fp2Trait>(
    generator: &Fp2,
    value: &Fp2,
    order: usize,
) -> (usize, u32) {
    let mut retval = FAILURE_RETVAL;

    let mut element = Fp2::ONE;
    let mut result = 0;
    for i in 0..order {
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
        .map(|exp| word_bn_to_byte_bn(&prime_power_to_bn_vartime(p, exp)))
        .collect::<Vec<_>>();

    let prime_order_subgroup_generator = generator.pow(
        &p_to_the_e_basis[e - 1].repr,
        p_to_the_e_basis[e - 1].bitlen,
    );
    assert_eq!(
        prime_order_subgroup_generator.equals(&Fp2::ONE),
        FAILURE_RETVAL,
        "g has order < {}^{}",
        p,
        e,
    );

    let mut partial_solutions = Vec::with_capacity(e);
    let mut partial_sum = vec![0];
    for i in 0..e {
        let partial_sum_bn = word_bn_to_byte_bn(&partial_sum);
        let r = *value
            * generator
                .pow(&partial_sum_bn.repr, partial_sum_bn.bitlen)
                .invert(); // TODO: can't we use the (much faster) conjugate since we're in a cyclotomic group?
        let u = r.pow(
            &p_to_the_e_basis[e - i - 1].repr,
            p_to_the_e_basis[e - i - 1].bitlen,
        );

        let (x, ok) = solve_dlp_small_prime_order(&prime_order_subgroup_generator, &u, p);
        partial_solutions.push(x);
        partial_sum = add_bn_vartime(
            &partial_sum,
            &mul_bn_by_u64_vartime(&byte_bn_to_word_bn(&p_to_the_e_basis[i]), x as u64),
        );
        retval &= ok;
    }

    (word_bn_to_byte_bn(&partial_sum), retval)
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

pub fn encrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    pub_key: &PubKey<Fp2>,
    message: &[u8],
) -> (Ciphertext<Fp2>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    // FIXME: where can I use vartime functions (i.e. operations on BigUint, gcd)? Where must things be constant-time?
    // TODO: add actual error handling

    let mut retval = SUCCESS_RETVAL;

    /* Sample scalars used for masking torsion points images or generating new kernels */

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = sample_random_element_mod(&pub_params.three_torsion_order);

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    // TODO: should this be full 2^a torsion, or effective 2^(a-2) torsion?
    let (omega, omega_inv) = sample_random_unit_mod(&pub_params.effective_two_torsion_order);

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    let D =
        sample_random_invertible_matrix_mod(&BigUint::from(5u8), &pub_params.five_torsion_order);

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>) and psi' (<R_A + [r] S_A>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r.repr,
        r.bitlen,
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.repr,
        r.bitlen,
    );

    // Apply sender's secret isogeny to 2^a-torsion basis to obtain their codomain curve E_B and basis image points (P_B, Q_B)
    let mut two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EB,
    );
    let two_torsion_basis_EB = BasisX::from_slice(&two_torsion_basis_EB);
    let (P_B, Q_B) = codomain_curve.lift_basis(&two_torsion_basis_EB);
    retval &= kernel_has_right_order;

    let masked_P_B = codomain_curve.mul(&P_B, &omega.repr, omega.bitlen);
    let masked_Q_B = codomain_curve.mul(&Q_B, &omega_inv.repr, omega_inv.bitlen);

    let masked_PQ_B = codomain_curve.sub(&masked_P_B, &masked_Q_B);

    let mut masked_two_torsion_basis_EB = [
        masked_P_B.to_pointx(),
        masked_Q_B.to_pointx(),
        masked_PQ_B.to_pointx(),
    ];
    PointX::batch_normalise(&mut masked_two_torsion_basis_EB);
    let masked_two_torsion_basis_EB = BasisX::from_slice(&masked_two_torsion_basis_EB);

    // Apply sender's secret isogeny to 5^c-torsion basis to obtain basis image points (X_B, Y_B)
    let mut five_torsion_basis_EB = pub_params.five_torsion_basis.to_array();
    let (codomain_curve_verif, kernel_has_right_order) =
        pub_params.starting_curve.three_isogeny_chain(
            &psi_kernel,
            pub_params.three_torsion_exp,
            &mut five_torsion_basis_EB,
        );
    let five_torsion_basis_EB = BasisX::from_slice(&five_torsion_basis_EB);
    let (X_B, Y_B) = codomain_curve_verif.lift_basis(&five_torsion_basis_EB);
    retval &= kernel_has_right_order;

    let masked_X_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, &D[0][0].repr, D[0][0].bitlen),
        &codomain_curve_verif.mul(&Y_B, &D[0][1].repr, D[0][1].bitlen),
    );
    let masked_Y_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, &D[1][0].repr, D[1][0].bitlen),
        &codomain_curve_verif.mul(&Y_B, &D[1][1].repr, D[1][1].bitlen),
    );


    let masked_XY_B = codomain_curve_verif.sub(&masked_X_B, &masked_Y_B);

    let mut masked_five_torsion_basis_EB = [
        masked_X_B.to_pointx(),
        masked_Y_B.to_pointx(),
        masked_XY_B.to_pointx(),
    ];
    PointX::batch_normalise(&mut masked_five_torsion_basis_EB);
    let masked_five_torsion_basis_EB = BasisX::from_slice(&masked_five_torsion_basis_EB);

    assert_eq!(
        codomain_curve
            .j_invariant()
            .equals(&codomain_curve_verif.j_invariant()),
        SUCCESS_RETVAL,
    );

    // Apply sender's secret parallel isogeny to receiver's masked 2^a-torsion basis image points to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB)
    let mut two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EAB,
    );
    let two_torsion_basis_EAB = BasisX::from_slice(&two_torsion_basis_EAB);
    let (P_AB, Q_AB) = shared_end_curve.lift_basis(&two_torsion_basis_EAB);
    retval &= kernel_has_right_order;

    let masked_P_AB = shared_end_curve.mul(&P_AB, &omega.repr, omega.bitlen);
    let masked_Q_AB = shared_end_curve.mul(&Q_AB, &omega_inv.repr, omega_inv.bitlen);

    let masked_PQ_AB = shared_end_curve.sub(&masked_P_AB, &masked_Q_AB);

    let mut masked_two_torsion_basis_EAB = [
        masked_P_AB.to_pointx(),
        masked_Q_AB.to_pointx(),
        masked_PQ_AB.to_pointx(),
    ];
    PointX::batch_normalise(&mut masked_two_torsion_basis_EAB);
    let masked_two_torsion_basis_EAB = BasisX::from_slice(&masked_two_torsion_basis_EAB);

    // Apply sender's secret parallel isogeny to receiver's masked 5^c-torsion basis image points to obtain shared secret (X_AB, Y_AB)
    let mut five_torsion_basis_EAB = pub_key.masked_five_torsion_basis_img.to_array();
    let (shared_end_curve_verif, kernel_has_right_order) =
        pub_key.codomain_curve.three_isogeny_chain(
            &psi_prime_kernel,
            pub_params.three_torsion_exp,
            &mut five_torsion_basis_EAB,
        );
    let five_torsion_basis_EAB = BasisX::from_slice(&five_torsion_basis_EAB);
    let (X_AB, Y_AB) = shared_end_curve_verif.lift_basis(&five_torsion_basis_EAB);
    retval &= kernel_has_right_order;

    let masked_X_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, &D[0][0].repr, D[0][0].bitlen),
        &shared_end_curve_verif.mul(&Y_AB, &D[0][1].repr, D[0][1].bitlen),
    );
    let masked_Y_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, &D[1][0].repr, D[1][0].bitlen),
        &shared_end_curve_verif.mul(&Y_AB, &D[1][1].repr, D[1][1].bitlen),
    );

    assert_eq!(
        shared_end_curve
            .j_invariant()
            .equals(&shared_end_curve_verif.j_invariant()),
        SUCCESS_RETVAL,
    );

    let mut kdf = Shake256::default();
    kdf.update(&masked_X_AB.to_pointx().x().encode());
    kdf.update(&masked_Y_AB.to_pointx().x().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut encrypted_message = vec![0u8; message.len()];
    let Ok(_) = one_time_pad.read(&mut encrypted_message) else {
        panic!("Could not read enough bytes from KDF");
    };
    for (encrypted_message_byte, message_byte) in encrypted_message.iter_mut().zip(message) {
        *encrypted_message_byte ^= message_byte;
    }

    let ct = Ciphertext {
        codomain_curve,
        masked_two_torsion_basis_EB,
        masked_five_torsion_basis_EB,
        shared_end_curve,
        masked_two_torsion_basis_EAB,
        encrypted_message,
    };

    (ct, retval)
}

pub fn decrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    prv_key: &PrvKey<Fp2>,
    ciphertext: &Ciphertext<Fp2>,
) -> (Vec<u8>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    // Factor that shows up in the application of the 2D-isogeny, from the dual that appears in the representation
    let dual_factor = &pub_params.effective_two_torsion_order - &prv_key.q;

    // Invert secret scalars, to neutralize their action on masked points we receive
    // TODO: should this be full 2^a torsion, or effective 2^(a-2) torsion?
    let alpha_inv = invert_element_mod(&prv_key.alpha, &pub_params.effective_two_torsion_order);
    let beta_inv = invert_element_mod(&prv_key.beta, &pub_params.effective_two_torsion_order);

    // Neutralize action of our own secret scalars on the masked 2^a-torsion basis for E_AB
    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);
    let unmasked_P_AB = ciphertext
        .shared_end_curve
        .mul(&P_AB, &alpha_inv.repr, alpha_inv.bitlen);
    let unmasked_Q_AB = ciphertext
        .shared_end_curve
        .mul(&Q_AB, &beta_inv.repr, beta_inv.bitlen);

    // Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let (P_B, Q_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EB);
    let mut deg_P_B = ciphertext
        .codomain_curve
        .mul(
            &P_B,
            &prv_key.q.to_bytes_le(),
            prv_key.q
                .bits()
                .try_into()
                .expect("Size in bits of the hidden degree q is too big to fit in a usize (we do not ever expect this to happen)"),
        );
    deg_P_B.set_neg();
    let mut deg_Q_B = ciphertext
        .codomain_curve
        .mul(
            &Q_B,
            &prv_key.q.to_bytes_le(),
            prv_key.q
                .bits()
                .try_into()
                .expect("Size in bits of the hidden degree q is too big to fit in a usize (we do not ever expect this to happen)"),
        );
    deg_Q_B.set_neg();

    let P1P2 = ProductPoint::new(&deg_P_B, &unmasked_P_AB);
    let Q1Q2 = ProductPoint::new(&deg_Q_B, &unmasked_Q_AB);

    // Compute Phi' on the masked 5^c-torsion for E_B
    // FIXME: requires points of order 2^(a+2)
    let domain = EllipticProduct::new(&ciphertext.codomain_curve, &ciphertext.shared_end_curve);
    let (X_B, Y_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_five_torsion_basis_EB);
    let XY_B = ciphertext.codomain_curve.sub(&X_B, &Y_B);
    let (aux_curves, five_torsion_basis_EB_on_aux_curve, ok) = domain.elliptic_product_isogeny(
        &P1P2,
        &Q1Q2,
        pub_params.effective_two_torsion_exp,
        &[
            ProductPoint::new(&X_B, &Point::INFINITY),
            ProductPoint::new(&Y_B, &Point::INFINITY),
            ProductPoint::new(&XY_B, &Point::INFINITY),
        ],
    );
    retval &= ok;

    // Generate random basis of the 5^c-torsion on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis(
        &ciphertext.shared_end_curve,
        &pub_params.five_torsion_order,
        &pub_params.five_torsion_cofactor,
    );
    let UV = ciphertext.shared_end_curve.sub(&U, &V);

    let (aux_curves_verif, five_torsion_basis_EAB_on_aux_curve, ok) = domain
        .elliptic_product_isogeny(
            &P1P2,
            &Q1Q2,
            pub_params.effective_two_torsion_exp,
            &[
                ProductPoint::new(&Point::INFINITY, &U),
                ProductPoint::new(&Point::INFINITY, &V),
                ProductPoint::new(&Point::INFINITY, &UV),
            ],
        );
    retval &= ok;

    assert_eq!(
        aux_curves
            .curves()
            .0
            .j_invariant()
            .equals(&aux_curves_verif.curves().0.j_invariant()),
        SUCCESS_RETVAL,
        "j-invariant of F1 in F1 x F2 does not match",
    );
    assert_eq!(
        aux_curves
            .curves()
            .1
            .j_invariant()
            .equals(&aux_curves_verif.curves().1.j_invariant()),
        SUCCESS_RETVAL,
        "j-invariant of F2 in F1 x F2 does not match",
    );

    /* Find change-of-basis matrix */

    // Compute pairs of point subtractions for later computing the pairings between them
    let mut X_intermediate_curve = five_torsion_basis_EB_on_aux_curve[0].points().0;
    let mut Y_intermediate_curve = five_torsion_basis_EB_on_aux_curve[1].points().0;
    let XY_intermediate_curve = five_torsion_basis_EB_on_aux_curve[2].points().0;
    Y_intermediate_curve.set_condneg(
        !aux_curves
            .curves()
            .0
            .sub(&X_intermediate_curve, &Y_intermediate_curve)
            .to_pointx()
            .equals(&XY_intermediate_curve.to_pointx()),
    );

    let mut U_intermediate_curve = five_torsion_basis_EAB_on_aux_curve[0].points().0;
    let mut V_intermediate_curve = five_torsion_basis_EAB_on_aux_curve[1].points().0;
    let UV_intermediate_curve = five_torsion_basis_EAB_on_aux_curve[2].points().0;
    V_intermediate_curve.set_condneg(
        !aux_curves_verif
            .curves()
            .0
            .sub(&U_intermediate_curve, &V_intermediate_curve)
            .to_pointx()
            .equals(&UV_intermediate_curve.to_pointx()),
    );

    let XV_intermediate_curve = aux_curves
        .curves()
        .0
        .sub(&X_intermediate_curve, &V_intermediate_curve);
    let XmU_intermediate_curve = aux_curves
        .curves()
        .0
        .add(&X_intermediate_curve, &U_intermediate_curve);

    let YV_intermediate_curve = aux_curves
        .curves()
        .0
        .sub(&Y_intermediate_curve, &V_intermediate_curve);
    let YmU_intermediate_curve = aux_curves
        .curves()
        .0
        .add(&Y_intermediate_curve, &U_intermediate_curve);

    // Compute the pairings e(U, V), e(X, V) = e(U, V)^x and e(X, -U) = e(U, V)^y,
    // e(Y, V) = e(U, V)^w and e(Y, -U) = e(U, V)^z
    // FIXME: Why does this direct way of computing the pairing not work?
    // let eUV_aux = aux_curves_verif.curves().0.weil_pairing(
    //     &U_intermediate_curve.to_pointx().x(),
    //     &V_intermediate_curve.to_pointx().x(),
    //     &UV_intermediate_curve.to_pointx().x(),
    //     &pub_params.five_torsion_order.to_bytes_le(),
    //     pub_params.five_torsion_order
    //         .bits()
    //         .try_into()
    //         .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    // );
    let eUV_aux = eUV_AB.pow(
        &dual_factor.to_bytes_le(),
        dual_factor
            .bits()
            .try_into()
            .expect("Size in bits of (2^a - q) is too big to fit into a usize (we don't expect this to ever happen)"),
    );

    // FIXME: none of the subsequent pairings are correct! This breaks everything!
    // I suspect a discrepancy between Sage's Weil pairing and the one here
    let eXV = aux_curves_verif.curves().0.weil_pairing(
        &X_intermediate_curve.to_pointx().x(),
        &V_intermediate_curve.to_pointx().x(),
        &XV_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let mU_intermediate_curve = -U_intermediate_curve;
    let eXmU = aux_curves_verif.curves().0.weil_pairing(
        &X_intermediate_curve.to_pointx().x(),
        &mU_intermediate_curve.to_pointx().x(),
        &XmU_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eYV = aux_curves_verif.curves().0.weil_pairing(
        &Y_intermediate_curve.to_pointx().x(),
        &V_intermediate_curve.to_pointx().x(),
        &YV_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eYmU = aux_curves_verif.curves().0.weil_pairing(
        &Y_intermediate_curve.to_pointx().x(),
        &mU_intermediate_curve.to_pointx().x(),
        &YmU_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    // Solve discrete logarithm between pairings to obtain expression of X' in terms of <U',V'>
    let (x, ok) = solve_dlp_small_prime_power_order(&eUV_aux, &eXV, 5, pub_params.five_torsion_exp);
    retval &= ok;
    let (y, ok) =
        solve_dlp_small_prime_power_order(&eUV_aux, &eXmU, 5, pub_params.five_torsion_exp);
    retval &= ok;
    let (w, ok) = solve_dlp_small_prime_power_order(&eUV_aux, &eYV, 5, pub_params.five_torsion_exp);
    retval &= ok;
    let (z, ok) =
        solve_dlp_small_prime_power_order(&eUV_aux, &eYmU, 5, pub_params.five_torsion_exp);
    retval &= ok;

    /* Decrypt message using one-time pad derived from shared secret */

    // Compute shared secret points (reusing temporary intermediate curve points as an optimization)
    ciphertext
        .shared_end_curve
        .mul_into(&mut X_intermediate_curve, &U, &x.repr, x.bitlen);
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_intermediate_curve, &V, &y.repr, y.bitlen);
    ciphertext.shared_end_curve.add_into(
        &mut U_intermediate_curve,
        &X_intermediate_curve,
        &Y_intermediate_curve,
    );
    ciphertext.shared_end_curve.mul_into(
        &mut V_intermediate_curve,
        &U_intermediate_curve,
        &dual_factor.to_bytes_le(),
        dual_factor
            .bits()
            .try_into()
            .expect("Size in bits of (2^a - q) is too big to fit in a usize (we do not ever expect this to happen)"),
    );
    let X_AB = ciphertext
        .shared_end_curve
        .mul(
            &V_intermediate_curve,
            &prv_key.delta.to_bytes_le(),
            prv_key.delta
                .bits()
                .try_into()
                .expect("Size in bits of delta is too big to fit in a usize (we do not ever expect this to happen)"),
        );

    ciphertext
        .shared_end_curve
        .mul_into(&mut X_intermediate_curve, &U, &w.repr, w.bitlen);
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_intermediate_curve, &V, &z.repr, z.bitlen);
    ciphertext.shared_end_curve.add_into(
        &mut U_intermediate_curve,
        &X_intermediate_curve,
        &Y_intermediate_curve,
    );
    ciphertext.shared_end_curve.mul_into(
        &mut V_intermediate_curve,
        &U_intermediate_curve,
        &dual_factor.to_bytes_le(),
        dual_factor
            .bits()
            .try_into()
            .expect("Size in bits of (2^a - q) is too big to fit in a usize (we do not ever expect this to happen)"),
    );
    let Y_AB = ciphertext
        .shared_end_curve
        .mul(
            &V_intermediate_curve,
            &prv_key.delta.to_bytes_le(),
            prv_key.delta
                .bits()
                .try_into()
                .expect("Size in bits of delta is too big to fit in a usize (we do not ever expect this to happen)"),
        );

    // Undo one-time pad of message
    let mut kdf = Shake256::default();
    kdf.update(&X_AB.to_pointx().x().encode());
    kdf.update(&Y_AB.to_pointx().x().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut message = vec![0u8; ciphertext.encrypted_message.len()];
    let Ok(_) = one_time_pad.read(&mut message) else {
        panic!("Could not read enough bytes from KDF");
    };
    println!("One-time pad used in decryption: {:?}", message);
    for (message_byte, encrypted_message_byte) in
        message.iter_mut().zip(&ciphertext.encrypted_message)
    {
        *message_byte ^= encrypted_message_byte;
    }

    (message, retval)
}

pub fn inke_encrypt<'a, Fp2: Fp2Trait>(
    pub_params: &InkePublicParams<Fp2>,
    pub_key: &InkePubKey<Fp2>,
    message: &[u8],
) -> (InkeCiphertext<Fp2>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    /* Sample scalars used for masking torsion points images or generating new kernels */

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = sample_random_element_mod(&pub_params.three_torsion_order);

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    let (omega, omega_inv) = sample_random_unit_mod(&pub_params.effective_two_torsion_order);

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>), psi' (<R_A + [r] S_A>) and psi'' (<R_A1 + [r] S_A1>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r.repr,
        r.bitlen,
    );
    let psi_prime_kernel = pub_key.intermediate_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img_intermediate,
        &r.repr,
        r.bitlen,
    );
    let psi_dblprime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.repr,
        r.bitlen,
    );

    // Apply sender's secret isogeny to 2^a-torsion basis to obtain their codomain curve E_B and basis image points (P_B, Q_B)
    let mut two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EB,
    );
    let two_torsion_basis_EB = BasisX::from_slice(&two_torsion_basis_EB);
    let (P_B, Q_B) = codomain_curve.lift_basis(&two_torsion_basis_EB);
    retval &= kernel_has_right_order;

    let masked_P_B = codomain_curve.mul(&P_B, &omega.repr, omega.bitlen);
    let masked_Q_B = codomain_curve.mul(&Q_B, &omega_inv.repr, omega_inv.bitlen);

    let masked_PQ_B = codomain_curve.sub(&masked_P_B, &masked_Q_B);

    let mut masked_two_torsion_basis_EB = [
        masked_P_B.to_pointx(),
        masked_Q_B.to_pointx(),
        masked_PQ_B.to_pointx(),
    ];
    PointX::batch_normalise(&mut masked_two_torsion_basis_EB);
    let masked_two_torsion_basis_EB = BasisX::from_slice(&masked_two_torsion_basis_EB);

    // Apply sender's secret parallel isogeny to receiver's masked 2^a-torsion basis image points to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB)
    let mut two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_dblprime_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EAB,
    );
    let two_torsion_basis_EAB = BasisX::from_slice(&two_torsion_basis_EAB);
    let (P_AB, Q_AB) = shared_end_curve.lift_basis(&two_torsion_basis_EAB);
    retval &= kernel_has_right_order;

    let masked_P_AB = shared_end_curve.mul(&P_AB, &omega.repr, omega.bitlen);
    let masked_Q_AB = shared_end_curve.mul(&Q_AB, &omega_inv.repr, omega_inv.bitlen);

    let masked_PQ_AB = shared_end_curve.sub(&masked_P_AB, &masked_Q_AB);

    let mut masked_two_torsion_basis_EAB = [
        masked_P_AB.to_pointx(),
        masked_Q_AB.to_pointx(),
        masked_PQ_AB.to_pointx(),
    ];
    PointX::batch_normalise(&mut masked_two_torsion_basis_EAB);
    let masked_two_torsion_basis_EAB = BasisX::from_slice(&masked_two_torsion_basis_EAB);

    // Compute codomain of sender's secret intermediate parallel isogeny to obtain shared secret curve
    let (secret_curve, _) = pub_key.intermediate_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut [],
    );

    // Compute shared secret from j-invariant of shared secret curve and encrypt message
    let mut kdf = Shake256::default();
    kdf.update(&secret_curve.j_invariant().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut encrypted_message = vec![0u8; message.len()];
    let Ok(_) = one_time_pad.read(&mut encrypted_message) else {
        panic!("Could not read enough bytes from KDF");
    };
    for (encrypted_message_byte, message_byte) in encrypted_message.iter_mut().zip(message) {
        *encrypted_message_byte ^= message_byte;
    }

    let ct = InkeCiphertext {
        codomain_curve,
        masked_two_torsion_basis_EB,
        shared_end_curve,
        masked_two_torsion_basis_EAB,
        encrypted_message,
    };

    (ct, retval)
}

pub fn inke_decrypt<Fp2: Fp2Trait>(
    pub_params: &InkePublicParams<Fp2>,
    prv_key: &InkePrvKey<Fp2>,
    ciphertext: &InkeCiphertext<Fp2>,
) -> (Vec<u8>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    // Invert secret scalars, to neutralize their action on masked points we receive
    // TODO: should this be full 2^a torsion, or effective 2^(a-2) torsion?
    let alpha_inv = invert_element_mod(&prv_key.alpha, &pub_params.effective_two_torsion_order);
    let beta_inv = invert_element_mod(&prv_key.beta, &pub_params.effective_two_torsion_order);

    // Neutralize action of our own secret scalars on the masked 2^a-torsion basis for E_AB
    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);
    let unmasked_P_AB = ciphertext
        .shared_end_curve
        .mul(&P_AB, &alpha_inv.repr, alpha_inv.bitlen);
    let unmasked_Q_AB = ciphertext
        .shared_end_curve
        .mul(&Q_AB, &beta_inv.repr, beta_inv.bitlen);

    // Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let (P_B, Q_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EB);
    let mut deg_P_B = ciphertext
        .codomain_curve
        .mul(
            &P_B,
            &prv_key.q.to_bytes_le(),
            prv_key.q
                .bits()
                .try_into()
                .expect("Size in bits of the hidden degree q is too big to fit in a usize (we do not ever expect this to happen)"),
        );
    deg_P_B.set_neg();
    let mut deg_Q_B = ciphertext
        .codomain_curve
        .mul(
            &Q_B,
            &prv_key.q.to_bytes_le(),
            prv_key.q
                .bits()
                .try_into()
                .expect("Size in bits of the hidden degree q is too big to fit in a usize (we do not ever expect this to happen)"),
        );
    deg_Q_B.set_neg();

    let P1P2 = ProductPoint::new(&deg_P_B, &unmasked_P_AB);
    let Q1Q2 = ProductPoint::new(&deg_Q_B, &unmasked_Q_AB);

    // Compute codomain curve pair of Phi', which contains the shared secret curve
    let domain = EllipticProduct::new(&ciphertext.codomain_curve, &ciphertext.shared_end_curve);
    let (aux_curves, _, ok) =
        domain.elliptic_product_isogeny(&P1P2, &Q1Q2, pub_params.effective_two_torsion_exp, &[]);
    retval &= ok;
    let secret_curve = aux_curves.curves().0;

    // Undo one-time pad of message
    let mut kdf = Shake256::default();
    kdf.update(&secret_curve.j_invariant().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut message = vec![0u8; ciphertext.encrypted_message.len()];
    let Ok(_) = one_time_pad.read(&mut message) else {
        panic!("Could not read enough bytes from KDF");
    };
    for (message_byte, encrypted_message_byte) in
        message.iter_mut().zip(&ciphertext.encrypted_message)
    {
        *message_byte ^= encrypted_message_byte;
    }

    (message, retval)
}
