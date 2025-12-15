use std::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};
use num_bigint::BigUint;
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _, XofReader as _},
};

use crate::{
    SUCCESS_RETVAL,
    bn::BigNum,
    dlp::solve_dlp_small_prime_power_order,
    rand::{
        sample_random_element_mod, sample_random_invertible_matrix_mod_prime_power,
        sample_random_torsion_basis_order_prime_power, sample_random_unit_mod_prime_power,
    },
};

pub struct PublicParams<Fp2: Fp2Trait> {
    pub starting_curve: Curve<Fp2>,
    pub full_two_torsion_order: BigNum,
    pub full_two_torsion_exp: usize,
    pub effective_two_torsion_order: BigNum,
    pub effective_two_torsion_exp: usize,
    pub three_torsion_order: BigNum,
    pub three_torsion_exp: usize,
    pub five_torsion_order: BigNum,
    pub five_torsion_exp: usize,
    pub five_torsion_cofactor: BigNum,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
    pub five_torsion_basis: BasisX<Fp2>,
    pub five_adic_basis: Vec<BigNum>,
}

// FIXME: represent scalars as their LE byte arrays and bitsize. Removes external dependency on num-bigint
pub struct PrvKey<Fp2: Fp2Trait> {
    pub q: BigNum,
    pub alpha: BigNum,
    pub beta: BigNum,
    pub delta: BigNum,
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

pub fn encrypt<Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    pub_key: &PubKey<Fp2>,
    message: &[u8],
) -> (Ciphertext<Fp2>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    // FIXME: where can I use vartime functions (i.e. operations on BigUint, gcd)? Where must things be constant-time?

    let mut retval = SUCCESS_RETVAL;

    /* Sample scalars used for masking torsion points images or generating new kernels */

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = sample_random_element_mod(&pub_params.three_torsion_order);

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    // TODO: should this be full 2^a torsion, or effective 2^(a-2) torsion?
    let omega1 = sample_random_unit_mod_prime_power(2, &pub_params.effective_two_torsion_order);
    let omega2 = sample_random_unit_mod_prime_power(2, &pub_params.effective_two_torsion_order);

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    let D = sample_random_invertible_matrix_mod_prime_power(5, &pub_params.five_torsion_order);

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>) and psi' (<R_A + [r] S_A>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        r.as_le_bytes(),
        r.nbits(),
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        r.as_le_bytes(),
        r.nbits(),
    );

    // Apply sender's secret isogeny to 2^a-torsion basis to obtain their codomain curve E_B and basis image points (P_B, Q_B)
    let mut two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EB,
    );
    // let [P_B, Q_B, ..] = &two_torsion_basis_EB;
    let two_torsion_basis_EB = BasisX::from_slice(&two_torsion_basis_EB);
    let (P_B, Q_B) = codomain_curve.lift_basis(&two_torsion_basis_EB);
    retval &= kernel_has_right_order;

    // let masked_xP_B = codomain_curve.xmul(P_B, &omega, omega_bitsize);
    // let masked_xQ_B = codomain_curve.xmul(Q_B, &omega_inv, omega_inv_bitsize);
    let masked_P_B = codomain_curve.mul(&P_B, omega1.as_le_bytes(), omega1.nbits());
    let masked_Q_B = codomain_curve.mul(&Q_B, omega2.as_le_bytes(), omega2.nbits());

    // let (masked_P_B, ok) = codomain_curve.lift_pointx(&masked_xP_B);
    // retval &= ok;
    // let (masked_Q_B, ok) = codomain_curve.lift_pointx(&masked_xQ_B);
    // retval &= ok;

    let masked_PQ_B = codomain_curve.sub(&masked_P_B, &masked_Q_B);

    // let masked_two_torsion_basis_EB =
    //     BasisX::from_points(&masked_xP_B, &masked_xQ_B, &masked_PQ_B.to_pointx());
    let masked_two_torsion_basis_EB = BasisX::from_points(
        &masked_P_B.to_pointx(),
        &masked_Q_B.to_pointx(),
        &masked_PQ_B.to_pointx(),
    );

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

    // let masked_xX_B = codomain_curve_verif.ladder_biscalar(
    //     &five_torsion_basis_EB,
    //     &D[(0, 0)],
    //     &D[(0, 1)],
    //     D_bitsize[(0, 0)],
    //     D_bitsize[(0, 1)],
    // );
    // let masked_xY_B = codomain_curve_verif.ladder_biscalar(
    //     &five_torsion_basis_EB,
    //     &D[(1, 0)],
    //     &D[(1, 1)],
    //     D_bitsize[(1, 0)],
    //     D_bitsize[(1, 1)],
    // );
    let masked_X_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, D[0][0].as_le_bytes(), D[0][0].nbits()),
        &codomain_curve_verif.mul(&Y_B, D[0][1].as_le_bytes(), D[0][1].nbits()),
    );
    let masked_Y_B = codomain_curve_verif.add(
        &codomain_curve_verif.mul(&X_B, D[1][0].as_le_bytes(), D[1][0].nbits()),
        &codomain_curve_verif.mul(&Y_B, D[1][1].as_le_bytes(), D[1][1].nbits()),
    );

    // let (masked_X_B, ok) = codomain_curve_verif.lift_pointx(&masked_xX_B);
    // retval &= ok;
    // let (masked_Y_B, ok) = codomain_curve_verif.lift_pointx(&masked_xY_B);
    // retval &= ok;

    let masked_XY_B = codomain_curve_verif.sub(&masked_X_B, &masked_Y_B);

    // let masked_five_torsion_basis_EB =
    //     BasisX::from_points(&masked_xX_B, &masked_xY_B, &masked_XY_B.to_pointx());
    let masked_five_torsion_basis_EB = BasisX::from_points(
        &masked_X_B.to_pointx(),
        &masked_Y_B.to_pointx(),
        &masked_XY_B.to_pointx(),
    );

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
    // let [xP_AB, xQ_AB, ..] = &two_torsion_basis_EAB;
    let two_torsion_basis_EAB = BasisX::from_slice(&two_torsion_basis_EAB);
    let (P_AB, Q_AB) = shared_end_curve.lift_basis(&two_torsion_basis_EAB);
    retval &= kernel_has_right_order;

    // let masked_xP_AB = shared_end_curve.xmul(xP_AB, &omega, omega_bitsize);
    // let masked_xQ_AB = shared_end_curve.xmul(xQ_AB, &omega_inv, omega_inv_bitsize);
    let masked_P_AB = shared_end_curve.mul(&P_AB, omega1.as_le_bytes(), omega1.nbits());
    let masked_Q_AB = shared_end_curve.mul(&Q_AB, omega2.as_le_bytes(), omega2.nbits());

    // let (masked_P_AB, ok) = shared_end_curve.lift_pointx(&masked_xP_AB);
    // retval &= ok;
    // let (masked_Q_AB, ok) = shared_end_curve.lift_pointx(&masked_xQ_AB);
    // retval &= ok;

    let masked_PQ_AB = shared_end_curve.sub(&masked_P_AB, &masked_Q_AB);

    // let masked_two_torsion_basis_EAB =
    //     BasisX::from_points(&masked_xP_AB, &masked_xQ_AB, &masked_PQ_AB.to_pointx());
    let masked_two_torsion_basis_EAB = BasisX::from_points(
        &masked_P_AB.to_pointx(),
        &masked_Q_AB.to_pointx(),
        &masked_PQ_AB.to_pointx(),
    );

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

    // let masked_xX_AB = shared_end_curve_verif.ladder_biscalar(
    //     &shared_secret,
    //     &D[(0, 0)],
    //     &D[(0, 1)],
    //     D_bitsize[(0, 0)],
    //     D_bitsize[(0, 1)],
    // );
    // let masked_xY_AB = shared_end_curve_verif.ladder_biscalar(
    //     &shared_secret,
    //     &D[(1, 0)],
    //     &D[(1, 1)],
    //     D_bitsize[(1, 0)],
    //     D_bitsize[(1, 1)],
    // );
    let masked_X_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, D[0][0].as_le_bytes(), D[0][0].nbits()),
        &shared_end_curve_verif.mul(&Y_AB, D[0][1].as_le_bytes(), D[0][1].nbits()),
    );
    let masked_Y_AB = shared_end_curve_verif.add(
        &shared_end_curve_verif.mul(&X_AB, D[1][0].as_le_bytes(), D[1][0].nbits()),
        &shared_end_curve_verif.mul(&Y_AB, D[1][1].as_le_bytes(), D[1][1].nbits()),
    );

    // let (masked_X_AB, ok) = shared_end_curve_verif.lift_pointx(&masked_xX_AB);
    // retval &= ok;
    // let (masked_Y_AB, ok) = shared_end_curve_verif.lift_pointx(&masked_xY_AB);
    // retval &= ok;

    // let masked_XY_AB = shared_end_curve_verif.sub(&masked_X_AB, &masked_Y_AB);

    // let shared_secret =
    //     BasisX::from_points(&masked_xX_AB, &masked_xY_AB, &masked_XY_AB.to_pointx());

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
    one_time_pad.read(&mut encrypted_message);
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
    let dual_factor = BigNum::new(
        &(BigUint::from_bytes_le(pub_params.effective_two_torsion_order.as_le_bytes())
            - BigUint::from_bytes_le(prv_key.q.as_le_bytes()))
        .to_u64_digits(),
    );

    // Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let (P_B, Q_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EB);

    let alpha_q = &prv_key.alpha * &prv_key.q;
    let mut scaled_P_B =
        ciphertext
            .codomain_curve
            .mul(&P_B, alpha_q.as_le_bytes(), alpha_q.nbits());
    scaled_P_B.set_neg();

    let beta_q = &prv_key.beta * &prv_key.q;
    let mut scaled_Q_B = ciphertext
        .codomain_curve
        .mul(&Q_B, beta_q.as_le_bytes(), beta_q.nbits());
    scaled_Q_B.set_neg();

    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);

    let P1P2 = ProductPoint::new(&scaled_P_B, &P_AB);
    let Q1Q2 = ProductPoint::new(&scaled_Q_B, &Q_AB);

    // Compute Phi' on the masked 5^c-torsion for E_B
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
    let (U, V, eUV_AB) = sample_random_torsion_basis_order_prime_power(
        &ciphertext.shared_end_curve,
        5,
        &pub_params.five_torsion_order,
        &pub_params.five_torsion_cofactor,
    );
    let UV = ciphertext.shared_end_curve.sub(&U, &V);

    // FIXME: do this only once, on an array of 4 image points ((X, 0), (Y, 0), (0, U), (0, V))
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
    let mut XY_aux_curve = five_torsion_basis_EB_on_aux_curve[2].points().0;
    Y_aux_curve.set_condneg(
        !aux_curve
            .sub(&X_aux_curve, &Y_aux_curve)
            .to_pointx()
            .equals(&XY_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = five_torsion_basis_EAB_on_aux_curve[0].points().0;
    let mut V_aux_curve = five_torsion_basis_EAB_on_aux_curve[1].points().0;
    let mut UV_aux_curve = five_torsion_basis_EAB_on_aux_curve[2].points().0;
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
    let eUV_aux = aux_curve.weil_pairing(
        &U_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &UV_aux_curve.to_pointx().x(),
        pub_params.five_torsion_order.as_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );
    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix that in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(prv_key.q.as_le_bytes(), prv_key.q.nbits());
    let eUV_aux_is_eUV_power_q = eUV_aux.equals(&eUV_power_q);

    // FIXME: none of the subsequent pairings are correct! This breaks everything!
    // I suspect a discrepancy between Sage's Weil pairing and the one here
    let eXV = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &XV_aux_curve.to_pointx().x(),
        pub_params.five_torsion_order.as_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );

    let mU_aux_curve = -U_aux_curve;
    let eXmU = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &mU_aux_curve.to_pointx().x(),
        &XmU_aux_curve.to_pointx().x(),
        pub_params.five_torsion_order.as_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );

    let eYV = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &YV_aux_curve.to_pointx().x(),
        pub_params.five_torsion_order.as_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );

    let eYmU = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &mU_aux_curve.to_pointx().x(),
        &YmU_aux_curve.to_pointx().x(),
        pub_params.five_torsion_order.as_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );

    // Solve discrete logarithm between pairings to obtain expression of X' in terms of <U',V'>
    let (x, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eXV,
        5,
        pub_params.five_torsion_exp,
        &pub_params.five_adic_basis,
    );
    retval &= ok;
    let (y, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eXmU,
        5,
        pub_params.five_torsion_exp,
        &pub_params.five_adic_basis,
    );
    retval &= ok;
    let (w, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eYV,
        5,
        pub_params.five_torsion_exp,
        &pub_params.five_adic_basis,
    );
    retval &= ok;
    let (z, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eYmU,
        5,
        pub_params.five_torsion_exp,
        &pub_params.five_adic_basis,
    );
    retval &= ok;

    /* Decrypt message using one-time pad derived from shared secret */

    // Compute shared secret points (reusing temporary intermediate curve points as an optimization)
    ciphertext
        .shared_end_curve
        .mul_into(&mut X_aux_curve, &U, x.as_le_bytes(), x.nbits());
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_aux_curve, &V, y.as_le_bytes(), y.nbits());
    ciphertext
        .shared_end_curve
        .add_into(&mut XY_aux_curve, &X_aux_curve, &Y_aux_curve);

    // [q] * X_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut U_aux_curve,
        &XY_aux_curve,
        prv_key.q.as_le_bytes(),
        prv_key.q.nbits(),
    );
    // [2^(a-2) - q] * X_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut V_aux_curve,
        &XY_aux_curve,
        dual_factor.as_le_bytes(),
        dual_factor.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    UV_aux_curve = U_aux_curve;
    UV_aux_curve.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

    let X_AB = ciphertext.shared_end_curve.mul(
        &UV_aux_curve,
        prv_key.delta.as_le_bytes(),
        prv_key.delta.nbits(),
    );

    ciphertext
        .shared_end_curve
        .mul_into(&mut X_aux_curve, &U, w.as_le_bytes(), w.nbits());
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_aux_curve, &V, z.as_le_bytes(), z.nbits());
    ciphertext
        .shared_end_curve
        .add_into(&mut XY_aux_curve, &X_aux_curve, &Y_aux_curve);

    // [q] * Y_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut U_aux_curve,
        &XY_aux_curve,
        prv_key.q.as_le_bytes(),
        prv_key.q.nbits(),
    );
    // [2^(a-2) - q] * Y_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut V_aux_curve,
        &XY_aux_curve,
        dual_factor.as_le_bytes(),
        dual_factor.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    UV_aux_curve = U_aux_curve;
    UV_aux_curve.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

    let Y_AB = ciphertext.shared_end_curve.mul(
        &UV_aux_curve,
        prv_key.delta.as_le_bytes(),
        prv_key.delta.nbits(),
    );

    // Undo one-time pad of message
    let mut kdf = Shake256::default();
    kdf.update(&X_AB.to_pointx().x().encode());
    kdf.update(&Y_AB.to_pointx().x().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut message = vec![0u8; ciphertext.encrypted_message.len()];
    one_time_pad.read(&mut message);
    for (message_byte, encrypted_message_byte) in
        message.iter_mut().zip(&ciphertext.encrypted_message)
    {
        *message_byte ^= encrypted_message_byte;
    }

    (message, retval)
}
