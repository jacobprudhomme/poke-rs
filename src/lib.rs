#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]

use std::{io::Read as _, marker::PhantomData};

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};
use num_bigint::BigUint;
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _},
};

use crate::{
    dlp::solve_dlp_small_prime_power_order,
    rand::{
        sample_random_element_mod, sample_random_invertible_matrix_mod,
        sample_random_torsion_basis, sample_random_unit_mod,
    },
    utilities::invert_element_mod,
};

mod bn;
pub mod dlp;
pub mod example_keypairs;
pub mod fields;
pub mod params;
pub mod rand;
mod utilities;

pub const SUCCESS_RETVAL: u32 = u32::MAX;
pub const FAILURE_RETVAL: u32 = u32::MIN;

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

// pub fn sample_2d_isogeny<Fp2: Fp2Trait>(
//     pub_params: &PublicParams<Fp2>,
// ) -> (
//     EllipticProduct<Fp2>,
//     ProductPoint<Fp2>,
//     ProductPoint<Fp2>,
//     u32,
// ) {
//     unimplemented!()
// }

// pub fn keygen<Fp2: Fp2Trait>(pub_params: &PublicParams<Fp2>) -> (PrvKey<Fp2>, PubKey<Fp2>, u32) {
//     let mut retval = SUCCESS_RETVAL;

//     /* Sample hidden degree of receiver's secret isogeny */
//     let mut rng = ndarray_rand::rand::thread_rng();
//     let ZERO = BigUint::from(0u8);
//     let ONE = BigUint::from(1u8);
//     let TWO = BigUint::from(2u8);
//     let THREE = BigUint::from(3u8);
//     let FIVE = BigUint::from(5u8);

//     let mut q = rng.gen_biguint_range(&ONE, &pub_params.effective_two_torsion_order);
//     while !(&q % &TWO != ZERO && &q % &THREE != ZERO && &q % &FIVE != ZERO) {
//         q = rng.gen_biguint_range(&ONE, &pub_params.effective_two_torsion_order);
//     }

//     /* Generate secret isogeny of degree q*(2^(a-2) - q) */
//     /* Sample scalars used for masking torsion points images */
//     let mut alpha = rng.gen_biguint_range(&ONE, &pub_params.full_two_torsion_order);
//     let mut alpha_inv = alpha.modinv(&pub_params.full_two_torsion_order);
//     while let None = alpha_inv {
//         // println!("alpha = {} is not invertible, retrying", alpha);
//         alpha = rng.gen_biguint_range(&ONE, &pub_params.full_two_torsion_order);
//         alpha_inv = alpha.modinv(&pub_params.full_two_torsion_order);
//     }
//     // println!();
//     let Some(alpha_inv) = alpha_inv else {
//         unreachable!();
//     };

//     let mut beta = rng.gen_biguint_range(&ONE, &pub_params.full_two_torsion_order);
//     let mut beta_inv = beta.modinv(&pub_params.full_two_torsion_order);
//     while let None = beta_inv {
//         // println!("beta = {} is not invertible, retrying", beta);
//         beta = rng.gen_biguint_range(&ONE, &pub_params.full_two_torsion_order);
//         beta_inv = beta.modinv(&pub_params.full_two_torsion_order);
//     }
//     // println!();
//     let Some(beta_inv) = beta_inv else {
//         unreachable!();
//     };

//     let mut gamma = rng.gen_biguint_range(&ONE, &pub_params.three_torsion_order);
//     let mut gamma_inv = gamma.modinv(&pub_params.three_torsion_order);
//     while let None = gamma_inv {
//         // println!("gamma = {} is not invertible, retrying", gamma);
//         gamma = rng.gen_biguint_range(&ONE, &pub_params.three_torsion_order);
//         gamma_inv = gamma.modinv(&pub_params.three_torsion_order);
//     }
//     // println!();
//     let Some(gamma_inv) = gamma_inv else {
//         unreachable!();
//     };

//     let mut delta = rng.gen_biguint_range(&ONE, &pub_params.five_torsion_order);
//     let mut delta_inv = delta.modinv(&pub_params.five_torsion_order);
//     while let None = delta_inv {
//         // println!("delta = {} is not invertible, retrying", delta);
//         delta = rng.gen_biguint_range(&ONE, &pub_params.five_torsion_order);
//         delta_inv = delta.modinv(&pub_params.five_torsion_order);
//     }
//     // println!();
//     let Some(delta_inv) = delta_inv else {
//         unreachable!();
//     };

//     /* Construct keypair */
//     // let prv_key = PrvKey {
//     //     q,
//     //     alpha,
//     //     beta,
//     //     delta,
//     //     _field: PhantomData,
//     // };

//     // let pub_key = PubKey {
//     //     codomain_curve: (),
//     //     masked_two_torsion_basis_img: (),
//     //     masked_three_torsion_basis_img: (),
//     //     masked_five_torsion_basis_img: (),
//     // };

//     // (prv_key, pub_key, retval)

//     unimplemented!()
// }

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
        &r.as_le_bytes(),
        r.nbits(),
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.as_le_bytes(),
        r.nbits(),
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

    let masked_P_B = codomain_curve.mul(&P_B, &omega.as_le_bytes(), omega.nbits());
    let masked_Q_B = codomain_curve.mul(&Q_B, &omega_inv.as_le_bytes(), omega_inv.nbits());

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
        &codomain_curve_verif.mul(&X_B, &D[0][0].as_le_bytes(), D[0][0].nbits()),
        &codomain_curve_verif.mul(&Y_B, &D[0][1].as_le_bytes(), D[0][1].nbits()),
    );
    let masked_Y_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, &D[1][0].as_le_bytes(), D[1][0].nbits()),
        &codomain_curve_verif.mul(&Y_B, &D[1][1].as_le_bytes(), D[1][1].nbits()),
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

    let masked_P_AB = shared_end_curve.mul(&P_AB, &omega.as_le_bytes(), omega.nbits());
    let masked_Q_AB = shared_end_curve.mul(&Q_AB, &omega_inv.as_le_bytes(), omega_inv.nbits());

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
        &shared_end_curve_verif.mul(&X_AB, &D[0][0].as_le_bytes(), D[0][0].nbits()),
        &shared_end_curve_verif.mul(&Y_AB, &D[0][1].as_le_bytes(), D[0][1].nbits()),
    );
    let masked_Y_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, &D[1][0].as_le_bytes(), D[1][0].nbits()),
        &shared_end_curve_verif.mul(&Y_AB, &D[1][1].as_le_bytes(), D[1][1].nbits()),
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

pub fn decrypt<Fp2: Fp2Trait>(
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
    let unmasked_P_AB =
        ciphertext
            .shared_end_curve
            .mul(&P_AB, &alpha_inv.as_le_bytes(), alpha_inv.nbits());
    let unmasked_Q_AB =
        ciphertext
            .shared_end_curve
            .mul(&Q_AB, &beta_inv.as_le_bytes(), beta_inv.nbits());

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
    let aux_curve = aux_curves.curves().0;

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
    let mut X_aux_curve = five_torsion_basis_EB_on_aux_curve[0].points().0;
    let mut Y_aux_curve = five_torsion_basis_EB_on_aux_curve[1].points().0;
    let XY_aux_curve = five_torsion_basis_EB_on_aux_curve[2].points().0;
    Y_aux_curve.set_condneg(
        !aux_curve
            .sub(&X_aux_curve, &Y_aux_curve)
            .to_pointx()
            .equals(&XY_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = five_torsion_basis_EAB_on_aux_curve[0].points().0;
    let mut V_aux_curve = five_torsion_basis_EAB_on_aux_curve[1].points().0;
    let UV_aux_curve = five_torsion_basis_EAB_on_aux_curve[2].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    let XV_aux_curve = aux_curve.sub(&X_aux_curve, &V_aux_curve);
    let XmU_aux_curve = aux_curve.add(&X_aux_curve, &U_aux_curve);

    let YV_aux_curve = aux_curve.sub(&Y_aux_curve, &V_aux_curve);
    let YmU_aux_curve = aux_curve.add(&Y_aux_curve, &U_aux_curve);

    // Compute the pairings e(U, V), e(X, V) = e(U, V)^x and e(X, -U) = e(U, V)^y,
    // e(Y, V) = e(U, V)^w and e(Y, -U) = e(U, V)^z
    // FIXME: Why does this direct way of computing the pairing not work?
    // let eUV_aux = aux_curve.weil_pairing(
    //     &U_aux_curve.to_pointx().x(),
    //     &V_aux_curve.to_pointx().x(),
    //     &UV_aux_curve.to_pointx().x(),
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
    let eXV = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &XV_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let mU_aux_curve = -U_aux_curve;
    let eXmU = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &mU_aux_curve.to_pointx().x(),
        &XmU_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eYV = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &YV_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_bytes_le(),
        pub_params.five_torsion_order
            .bits()
            .try_into()
            .expect("Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)"),
    );

    let eYmU = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &mU_aux_curve.to_pointx().x(),
        &YmU_aux_curve.to_pointx().x(),
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
        .mul_into(&mut X_aux_curve, &U, &x.as_le_bytes(), x.nbits());
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_aux_curve, &V, &y.as_le_bytes(), y.nbits());
    ciphertext
        .shared_end_curve
        .add_into(&mut U_aux_curve, &X_aux_curve, &Y_aux_curve);
    ciphertext.shared_end_curve.mul_into(
        &mut V_aux_curve,
        &U_aux_curve,
        &dual_factor.to_bytes_le(),
        dual_factor
            .bits()
            .try_into()
            .expect("Size in bits of (2^a - q) is too big to fit in a usize (we do not ever expect this to happen)"),
    );
    let X_AB = ciphertext
        .shared_end_curve
        .mul(
            &V_aux_curve,
            &prv_key.delta.to_bytes_le(),
            prv_key.delta
                .bits()
                .try_into()
                .expect("Size in bits of delta is too big to fit in a usize (we do not ever expect this to happen)"),
        );

    ciphertext
        .shared_end_curve
        .mul_into(&mut X_aux_curve, &U, &w.as_le_bytes(), w.nbits());
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_aux_curve, &V, &z.as_le_bytes(), z.nbits());
    ciphertext
        .shared_end_curve
        .add_into(&mut U_aux_curve, &X_aux_curve, &Y_aux_curve);
    ciphertext.shared_end_curve.mul_into(
        &mut V_aux_curve,
        &U_aux_curve,
        &dual_factor.to_bytes_le(),
        dual_factor
            .bits()
            .try_into()
            .expect("Size in bits of (2^a - q) is too big to fit in a usize (we do not ever expect this to happen)"),
    );
    let Y_AB = ciphertext
        .shared_end_curve
        .mul(
            &V_aux_curve,
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
        &r.as_le_bytes(),
        r.nbits(),
    );
    let psi_prime_kernel = pub_key.intermediate_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img_intermediate,
        &r.as_le_bytes(),
        r.nbits(),
    );
    let psi_dblprime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.as_le_bytes(),
        r.nbits(),
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

    let masked_P_B = codomain_curve.mul(&P_B, &omega.as_le_bytes(), omega.nbits());
    let masked_Q_B = codomain_curve.mul(&Q_B, &omega_inv.as_le_bytes(), omega_inv.nbits());

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

    let masked_P_AB = shared_end_curve.mul(&P_AB, &omega.as_le_bytes(), omega.nbits());
    let masked_Q_AB = shared_end_curve.mul(&Q_AB, &omega_inv.as_le_bytes(), omega_inv.nbits());

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
    let unmasked_P_AB =
        ciphertext
            .shared_end_curve
            .mul(&P_AB, &alpha_inv.as_le_bytes(), alpha_inv.nbits());
    let unmasked_Q_AB =
        ciphertext
            .shared_end_curve
            .mul(&Q_AB, &beta_inv.as_le_bytes(), beta_inv.nbits());

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
