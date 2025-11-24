#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]

use std::{io::Read as _, marker::PhantomData};

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
    utilities::bn::{
        add_bn_vartime, bn_bit_length_vartime, mul_bn_by_u64_vartime, prime_power_to_bn_vartime,
    },
};
use ndarray::Array2;
use ndarray_rand::{RandomExt as _, rand::distributions::Uniform};
use num_bigint::{BigUint, RandBigInt as _};
use num_integer::gcd;
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

// FIXME: represent scalars as their LE byte arrays and bitsize. Removes external dependency on num-bigint
pub struct PrvKey<Fp2: Fp2Trait> {
    pub q: BigUint,
    pub alpha: BigUint,
    pub beta: BigUint,
    pub delta: BigUint,
    pub _field: PhantomData<Fp2>,
}

pub struct PubKey<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_img: BasisX<Fp2>,
    pub masked_three_torsion_basis_img: BasisX<Fp2>,
    pub masked_five_torsion_basis_img: BasisX<Fp2>,
}

pub struct Ciphertext<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EB: BasisX<Fp2>,
    pub masked_five_torsion_basis_EB: BasisX<Fp2>,
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

pub fn solve_dlp_small_prime_order<Fp2: Fp2Trait>(
    generator: &Fp2,
    value: &Fp2,
    order: usize,
) -> (usize, u32) {
    let mut retval = FAILURE_RETVAL;

    println!("\nSolving in subgroup of order {}", order);
    println!("({}) = ({})^x", value, generator);

    let mut element = Fp2::ONE;
    let mut result = 0;
    for i in 0..order {
        let found_log = value.equals(&element);
        println!(
            "i = {}, {}",
            i,
            if found_log == SUCCESS_RETVAL {
                "FOUND"
            } else {
                "NOT FOUND"
            }
        );
        result |= i & (((found_log as usize) << 32) | found_log as usize);
        retval |= found_log;

        element *= *generator;
    }

    println!("Result: {}", result);

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

    println!("\nSolving in subgroup of order {}^{}", p, e);
    println!("({}) = ({})^x", value, generator);

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

    println!(
        "x = {} in ({}) = ({})^x",
        BigUint::from_bytes_le(&word_bn_to_byte_bn(&partial_sum).repr),
        value,
        generator,
    );

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
            println!(
                "Generated point U is in the torsion subgroup of the other cofactors of p + 1. Trying a new point"
            );
            continue;
        }

        let U_saturated = curve.mul(
            &U_in_torsion_subgroup,
            &reduced_torsion_subgroup_order,
            reduced_torsion_subgroup_order_bitsize,
        );
        if U_saturated.is_zero() == SUCCESS_RETVAL {
            println!("Generated point U has order < (pi)^(ei). Trying a new point");
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
            println!(
                "Generated point V is in the torsion subgroup of the other cofactors. Trying a new point"
            );
            continue;
        }

        let V_saturated = curve.mul(
            &V_in_torsion_subgroup,
            &reduced_torsion_subgroup_order,
            reduced_torsion_subgroup_order_bitsize,
        );
        if V_saturated.is_zero() == SUCCESS_RETVAL {
            println!("Generated point V has order < (pi)^(ei). Trying a new point");
            continue;
        }

        let V = V_in_torsion_subgroup;
        let UV = curve.sub(&U, &V);

        /* Check point is linearly independent to U */

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
            println!("e(U, V) does not have full multiplicative order");
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
            println!("e(U, V) has multiplicative order > (pi)^(ei)");
            continue;
        }

        break (V, eUV);
    };

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

    let mut rng = ndarray_rand::rand::thread_rng();
    let ONE = BigUint::from(1u8);

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = rng.gen_biguint_below(&pub_params.three_torsion_order); // FIXME: what happens if this is 0?
    let r_bitsize =
        r.bits().try_into().expect("Size in bits of the scalar r is too big to fit in a usize (we do not ever expect this to happen)");
    let r = r.to_bytes_le();

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    let mut omega = rng.gen_biguint_range(&ONE, &pub_params.effective_two_torsion_order);
    let mut omega_inv = omega.modinv(&pub_params.effective_two_torsion_order);
    while let None = omega_inv {
        println!("omega = {} is not invertible, retrying", omega);
        omega = rng.gen_biguint_range(&ONE, &pub_params.effective_two_torsion_order);
        omega_inv = omega.modinv(&pub_params.effective_two_torsion_order);
    }
    println!();
    let Some(omega_inv) = omega_inv else {
        unreachable!();
    };
    let omega_bitsize =
        omega.bits().try_into().expect("Size in bits of the scalar omega is too big to fit in a usize (we do not ever expect this to happen)");
    let omega_inv_bitsize =
        omega_inv.bits().try_into().expect("Size in bits of the scalar 1/omega is too big to fit in a usize (we do not ever expect this to happen)");
    let omega = omega.to_bytes_le();
    let omega_inv = omega_inv.to_bytes_le();

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    // FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
    let mut D = Array2::random_using(
        (2, 2),
        Uniform::new(BigUint::ZERO, &pub_params.five_torsion_order),
        &mut rng,
    );
    let mut det = (((&D[(0, 0)] * &D[(1, 1)]) % &pub_params.five_torsion_order)
        + (&pub_params.five_torsion_order
            - ((&D[(0, 1)] * &D[(1, 0)]) % &pub_params.five_torsion_order)))
        % &pub_params.five_torsion_order;
    let mut det_gcd = gcd(det.clone(), pub_params.five_torsion_order.clone()); // TODO: look into a borrowing GCD function
    while det_gcd != ONE {
        println!("det(D) = {}, gcd(det(D), 5^c) = {}, retrying", det, det_gcd);
        D = Array2::random_using(
            (2, 2),
            Uniform::new(BigUint::ZERO, &pub_params.five_torsion_order),
            &mut rng,
        );
        det = (((&D[(0, 0)] * &D[(1, 1)]) % &pub_params.five_torsion_order)
            + (&pub_params.five_torsion_order
                - ((&D[(0, 1)] * &D[(1, 0)]) % &pub_params.five_torsion_order)))
            % &pub_params.five_torsion_order;
        det_gcd = gcd(det.clone(), pub_params.five_torsion_order.clone()); // TODO: look into a borrowing GCD function
    }
    println!();
    let D_bitsize = D.map(|elem| {
        TryInto::<usize>::try_into(elem.bits())
            .expect("Size in bits of the scalar is too big to fit in a usize (we do not ever expect this to happen)")
    });
    let D = D.map(|elem| elem.to_bytes_le());

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>) and psi' (<R_A + [r] S_A>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r,
        r_bitsize,
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r,
        r_bitsize,
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
    println!(
        "Successful execution after applying psi to 2^a-torsion: {}",
        retval == SUCCESS_RETVAL,
    );

    let masked_P_B = codomain_curve.mul(&P_B, &omega, omega_bitsize);
    let masked_Q_B = codomain_curve.mul(&Q_B, &omega_inv, omega_inv_bitsize);

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
    println!(
        "Successful execution after applying psi to 5^c-torsion: {}",
        retval == SUCCESS_RETVAL,
    );

    let masked_X_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, &D[(0, 0)], D_bitsize[(0, 0)]),
        &codomain_curve_verif.mul(&Y_B, &D[(0, 1)], D_bitsize[(0, 1)]),
    );
    let masked_Y_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, &D[(1, 0)], D_bitsize[(1, 0)]),
        &codomain_curve_verif.mul(&Y_B, &D[(1, 1)], D_bitsize[(1, 1)]),
    );


    let masked_XY_B = codomain_curve_verif.sub(&masked_X_B, &masked_Y_B);

    let mut masked_five_torsion_basis_EB = [
        masked_X_B.to_pointx(),
        masked_Y_B.to_pointx(),
        masked_XY_B.to_pointx(),
    ];
    PointX::batch_normalise(&mut masked_five_torsion_basis_EB);
    let masked_five_torsion_basis_EB = BasisX::from_slice(&masked_five_torsion_basis_EB);

    println!("j-invariant for sender's codomain curve:");
    println!("{}", codomain_curve.j_invariant());
    println!("{}\n", codomain_curve_verif.j_invariant());
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
    println!(
        "Successful execution after applying psi' to 2^a-torsion: {}",
        retval == SUCCESS_RETVAL,
    );

    let masked_P_AB = shared_end_curve.mul(&P_AB, &omega, omega_bitsize);
    let masked_Q_AB = shared_end_curve.mul(&Q_AB, &omega_inv, omega_inv_bitsize);

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
    println!(
        "Successful execution after applying psi' to 5^c-torsion: {}",
        retval == SUCCESS_RETVAL,
    );

    let masked_X_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, &D[(0, 0)], D_bitsize[(0, 0)]),
        &shared_end_curve_verif.mul(&Y_AB, &D[(0, 1)], D_bitsize[(0, 1)]),
    );
    let masked_Y_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, &D[(1, 0)], D_bitsize[(1, 0)]),
        &shared_end_curve_verif.mul(&Y_AB, &D[(1, 1)], D_bitsize[(1, 1)]),
    );

    println!("j-invariant for the shared end curve:");
    println!("{}", shared_end_curve.j_invariant());
    println!("{}\n", shared_end_curve_verif.j_invariant());
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


    // Invert secret scalars, to neutralize their action on masked points we receive
    let alpha_inv = prv_key.alpha.modinv(&pub_params.full_two_torsion_order);
    let Some(alpha_inv) = alpha_inv else {
        unreachable!();
    };
    let alpha_inv_bitsize =
        alpha_inv.bits().try_into().expect("Size in bits of the scalar 1/alpha is too big to fit in a usize (we do not ever expect this to happen)");
    let alpha_inv = alpha_inv.to_bytes_le();

    let beta_inv = prv_key.beta.modinv(&pub_params.full_two_torsion_order);
    let Some(beta_inv) = beta_inv else {
        unreachable!();
    };
    let beta_inv_bitsize =
        beta_inv.bits().try_into().expect("Size in bits of the scalar 1/beta is too big to fit in a usize (we do not ever expect this to happen)");
    let beta_inv = beta_inv.to_bytes_le();

    // Neutralize action of our own secret scalars on the masked 2^a-torsion basis for E_AB
    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);
    let unmasked_P_AB = ciphertext
        .shared_end_curve
        .mul(&P_AB, &alpha_inv, alpha_inv_bitsize);
    let unmasked_Q_AB = ciphertext
        .shared_end_curve
        .mul(&Q_AB, &beta_inv, beta_inv_bitsize);

    // Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let (P_B, Q_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EB);
    let generator_point1_B = ciphertext
        .codomain_curve
        .mul(
            &P_B,
            &prv_key.q.to_bytes_le(),
            prv_key.q
                .bits()
                .try_into()
                .expect("Size in bits of the hidden degree q is too big to fit in a usize (we do not ever expect this to happen)"),
        );
    let generator_point2_B = ciphertext
        .codomain_curve
        .mul(
            &Q_B,
            &prv_key.q.to_bytes_le(),
            prv_key.q
                .bits()
                .try_into()
                .expect("Size in bits of the hidden degree q is too big to fit in a usize (we do not ever expect this to happen)"),
        );

    let kernel_generator_point1 = ProductPoint::new(&generator_point1_B, &unmasked_P_AB);
    let kernel_generator_point2 = ProductPoint::new(&generator_point2_B, &unmasked_Q_AB);

    // Compute Phi' on the masked 5^c-torsion for E_B
    // FIXME: requires points of order 2^(a+2)
    let domain = EllipticProduct::new(&ciphertext.codomain_curve, &ciphertext.shared_end_curve);
    let (X_B, Y_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_five_torsion_basis_EB);
    let (intermediate_curves, five_torsion_basis_intermediate_curves_left, ok) = domain
        .elliptic_product_isogeny(
            &kernel_generator_point1,
            &kernel_generator_point2,
            pub_params.effective_two_torsion_exp,
            &[
                ProductPoint::new(&X_B, &Point::INFINITY),
                ProductPoint::new(&Y_B, &Point::INFINITY),
            ],
        );
    retval &= ok;
    println!(
        "Successful execution after applying Phi' to 5^c-torsion on E_B: {}",
        retval == SUCCESS_RETVAL,
    );

    // Generate random basis of the 5^c-torsion on E_AB
    let (U, V, _) = sample_random_torsion_basis(
        &ciphertext.shared_end_curve,
        &pub_params.five_torsion_order,
        &pub_params.five_torsion_cofactor,
    );

    let (intermediate_curves_verif, five_torsion_basis_intermediate_curve_right, ok) = domain
        .elliptic_product_isogeny(
            &kernel_generator_point1,
            &kernel_generator_point2,
            pub_params.effective_two_torsion_exp,
            &[
                ProductPoint::new(&Point::INFINITY, &U),
                ProductPoint::new(&Point::INFINITY, &V),
            ],
        );
    retval &= ok;

    /* Find change-of-basis matrix */

    // Compute pairs of point subtractions for later computing the pairings between them
    let mut X_intermediate_curve = five_torsion_basis_intermediate_curves_left[0].points().0;
    let mut Y_intermediate_curve = five_torsion_basis_intermediate_curves_left[1].points().0;

    let mut U_intermediate_curve = five_torsion_basis_intermediate_curve_right[0].points().0;
    let mut V_intermediate_curve = five_torsion_basis_intermediate_curve_right[1].points().0;
    let UV_intermediate_curve = intermediate_curves_verif
        .curves()
        .0
        .sub(&U_intermediate_curve, &V_intermediate_curve);

    let XV_intermediate_curve = intermediate_curves
        .curves()
        .0
        .sub(&X_intermediate_curve, &V_intermediate_curve);
    let XmU_intermediate_curve = intermediate_curves
        .curves()
        .0
        .add(&X_intermediate_curve, &U_intermediate_curve);

    let YV_intermediate_curve = intermediate_curves
        .curves()
        .0
        .sub(&Y_intermediate_curve, &V_intermediate_curve);
    let YmU_intermediate_curve = intermediate_curves
        .curves()
        .0
        .add(&Y_intermediate_curve, &U_intermediate_curve);

    // Compute the pairings e(U, V), e(X, V) = e(U, V)^x and e(X, -U) = e(U, V)^y,
    // e(Y, V) = e(U, V)^w and e(Y, -U) = e(U, V)^z
    let eUV = intermediate_curves_verif.curves().0.weil_pairing(
        &U_intermediate_curve.to_pointx().x(),
        &V_intermediate_curve.to_pointx().x(),
        &UV_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eXV = intermediate_curves_verif.curves().0.weil_pairing(
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
    let eXmU = intermediate_curves_verif.curves().0.weil_pairing(
        &X_intermediate_curve.to_pointx().x(),
        &mU_intermediate_curve.to_pointx().x(),
        &XmU_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eYV = intermediate_curves_verif.curves().0.weil_pairing(
        &Y_intermediate_curve.to_pointx().x(),
        &V_intermediate_curve.to_pointx().x(),
        &YV_intermediate_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eYmU = intermediate_curves_verif.curves().0.weil_pairing(
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
    let (x, ok) = solve_dlp_small_prime_power_order(&eUV, &eXV, 5, pub_params.five_torsion_exp);
    retval &= ok;
    let (y, ok) = solve_dlp_small_prime_power_order(&eUV, &eXmU, 5, pub_params.five_torsion_exp);
    retval &= ok;
    let (w, ok) = solve_dlp_small_prime_power_order(&eUV, &eYV, 5, pub_params.five_torsion_exp);
    retval &= ok;
    let (z, ok) = solve_dlp_small_prime_power_order(&eUV, &eYmU, 5, pub_params.five_torsion_exp);
    retval &= ok;

    /* Decrypt message using one-time pad derived from shared secret */

    // Compute shared secret points (reusing temporary intermediate curve points as an optimization)
    let undoing_factor = &pub_params.effective_two_torsion_order - &prv_key.q;

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
        &undoing_factor.to_bytes_le(),
        undoing_factor
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
        &undoing_factor.to_bytes_le(),
        undoing_factor
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
    for (message_byte, encrypted_message_byte) in
        message.iter_mut().zip(&ciphertext.encrypted_message)
    {
        *message_byte ^= encrypted_message_byte;
    }

    (message, retval)
}
