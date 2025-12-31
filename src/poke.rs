use std::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _, XofReader as _},
};

use crate::{
    SUCCESS_RETVAL,
    bn::BigNum,
    dlp::solve_dlp_small_prime_power_order,
    masking::{
        mask_basisx_by_diagonal_scalars, mask_basisx_by_diagonal_scalars_points_only,
        mask_basisx_by_scalar_matrix, mask_basisx_by_scalar_matrix_pointx_only,
    },
    rand::{
        sample_random_element_mod, sample_random_invertible_matrix_mod_prime_power,
        sample_random_torsion_basis_order_prime_power, sample_random_unit_mod_prime_power,
    },
};

pub struct PublicParams<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_COF: usize,
    const NUM_WORDS_P: usize,
> {
    pub starting_curve: Curve<Fp2>,
    pub full_two_torsion_order: BigNum<NUM_WORDS_2>,
    pub full_two_torsion_exp: usize,
    pub effective_two_torsion_order: BigNum<NUM_WORDS_2>,
    pub effective_two_torsion_exp: usize,
    pub three_torsion_order: BigNum<NUM_WORDS_3>,
    pub three_torsion_exp: usize,
    pub five_torsion_order: BigNum<NUM_WORDS_5>,
    pub five_torsion_exp: usize,
    pub five_torsion_cofactor: BigNum<NUM_WORDS_COF>,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
    pub five_torsion_basis: BasisX<Fp2>,
    pub five_adic_basis: Vec<BigNum<NUM_WORDS_5>>,
}

pub struct PrvKey<Fp2: Fp2Trait, const NUM_WORDS_2: usize, const NUM_WORDS_5: usize> {
    pub q: BigNum<NUM_WORDS_2>,
    pub alpha: BigNum<NUM_WORDS_2>,
    pub beta: BigNum<NUM_WORDS_2>,
    pub delta: BigNum<NUM_WORDS_5>,
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

pub fn encrypt<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_COF: usize,
    const NUM_WORDS_P: usize,
>(
    pub_params: &PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_COF,
        NUM_WORDS_P,
    >,
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
    // FIXME: Should this be full 2^a torsion, or effective 2^(a-2) torsion?
    let omega1 = sample_random_unit_mod_prime_power(2, &pub_params.effective_two_torsion_order);
    let omega2 = sample_random_unit_mod_prime_power(2, &pub_params.effective_two_torsion_order);

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    let D = sample_random_invertible_matrix_mod_prime_power(5, &pub_params.five_torsion_order);

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>) and psi' (<R_A + [r] S_A>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r.to_le_bytes(),
        r.nbits(),
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.to_le_bytes(),
        r.nbits(),
    );

    // Apply sender's secret isogeny to 2^a- and 5^c-torsion bases to obtain their codomain curve E_B and basis image points (P_B, Q_B), (X_B, Y_B)
    let two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let five_torsion_basis_EB = pub_params.five_torsion_basis.to_array();
    let mut torsion_bases_EB = [PointX::INFINITY; 6];
    torsion_bases_EB[..3].copy_from_slice(&two_torsion_basis_EB);
    torsion_bases_EB[3..].copy_from_slice(&five_torsion_basis_EB);

    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut torsion_bases_EB,
    );
    retval &= kernel_has_right_order;

    let two_torsion_basis_EB = BasisX::from_slice(&torsion_bases_EB[..3]);
    let masked_two_torsion_basis_EB =
        mask_basisx_by_diagonal_scalars(&codomain_curve, &two_torsion_basis_EB, &omega1, &omega2);

    let five_torsion_basis_EB = BasisX::from_slice(&torsion_bases_EB[3..]);
    let masked_five_torsion_basis_EB =
        mask_basisx_by_scalar_matrix(&codomain_curve, &five_torsion_basis_EB, &D);

    // Apply sender's secret parallel isogeny to receiver's masked 2^a- and 5^c-torsion basis image points
    // to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB), (X_AB, Y_AB)
    let two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let five_torsion_basis_EAB = pub_key.masked_five_torsion_basis_img.to_array();
    let mut torsion_bases_EAB = [PointX::INFINITY; 6];
    torsion_bases_EAB[..3].copy_from_slice(&two_torsion_basis_EAB);
    torsion_bases_EAB[3..].copy_from_slice(&five_torsion_basis_EAB);

    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut torsion_bases_EAB,
    );
    retval &= kernel_has_right_order;

    let two_torsion_basis_EAB = BasisX::from_slice(&torsion_bases_EAB[..3]);
    let masked_two_torsion_basis_EAB = mask_basisx_by_diagonal_scalars(
        &shared_end_curve,
        &two_torsion_basis_EAB,
        &omega1,
        &omega2,
    );

    let five_torsion_basis_EAB = BasisX::from_slice(&torsion_bases_EAB[3..]);
    let (masked_xX_AB, masked_xY_AB) =
        mask_basisx_by_scalar_matrix_pointx_only(&shared_end_curve, &five_torsion_basis_EAB, &D);

    let mut kdf = Shake256::default();
    kdf.update(&masked_xX_AB.x().encode());
    kdf.update(&masked_xY_AB.x().encode());
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

pub fn decrypt<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_COF: usize,
    const NUM_WORDS_P: usize,
>(
    pub_params: &PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_COF,
        NUM_WORDS_P,
    >,
    prv_key: &PrvKey<Fp2, NUM_WORDS_2, NUM_WORDS_5>,
    ciphertext: &Ciphertext<Fp2>,
) -> (Vec<u8>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    // Factor that shows up in the application of the 2D-isogeny, from the dual that appears in the representation
    let dual_factor = &pub_params.effective_two_torsion_order - &prv_key.q;

    // Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let alpha_q = prv_key
        .alpha
        .mul_mod_power_of_two(&prv_key.q, &pub_params.full_two_torsion_order);
    let beta_q = prv_key
        .beta
        .mul_mod_power_of_two(&prv_key.q, &pub_params.full_two_torsion_order);

    let (scaled_P_B, scaled_Q_B) = mask_basisx_by_diagonal_scalars_points_only(
        &ciphertext.codomain_curve,
        &ciphertext.masked_two_torsion_basis_EB,
        &alpha_q,
        &beta_q,
    );

    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);

    let domain = EllipticProduct::new(&ciphertext.codomain_curve, &ciphertext.shared_end_curve);
    let P1P2 = ProductPoint::new(&scaled_P_B, &P_AB);
    let Q1Q2 = ProductPoint::new(&scaled_Q_B, &Q_AB);

    let (X_B, Y_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_five_torsion_basis_EB);
    let XY_B = ciphertext.codomain_curve.sub(&X_B, &Y_B);

    // Generate random basis of the 5^c-torsion on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis_order_prime_power(
        &ciphertext.shared_end_curve,
        5,
        &pub_params.five_torsion_order,
        &pub_params.five_torsion_cofactor,
    );
    let UV = ciphertext.shared_end_curve.sub(&U, &V);

    // Compute Phi' on the masked 5^c-torsion for E_B and a random 5^c-torsion basis for E_AB
    let (aux_curves, torsion_bases_aux_curves, ok) = domain.elliptic_product_isogeny(
        &P1P2,
        &Q1Q2,
        pub_params.effective_two_torsion_exp,
        &[
            ProductPoint::new(&X_B, &Point::INFINITY),
            ProductPoint::new(&Y_B, &Point::INFINITY),
            ProductPoint::new(&XY_B, &Point::INFINITY),
            ProductPoint::new(&Point::INFINITY, &U),
            ProductPoint::new(&Point::INFINITY, &V),
            ProductPoint::new(&Point::INFINITY, &UV),
        ],
    );
    retval &= ok;
    let aux_curve = aux_curves.curves().0;

    /* Find change-of-basis matrix */

    // Correct the pairs of image points to overall sign
    let mut X_aux_curve = torsion_bases_aux_curves[0].points().0;
    let mut Y_aux_curve = torsion_bases_aux_curves[1].points().0;
    let mut XY_aux_curve = torsion_bases_aux_curves[2].points().0;
    Y_aux_curve.set_condneg(
        !aux_curve
            .sub(&X_aux_curve, &Y_aux_curve)
            .to_pointx()
            .equals(&XY_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = torsion_bases_aux_curves[3].points().0;
    let mut V_aux_curve = torsion_bases_aux_curves[4].points().0;
    let mut UV_aux_curve = torsion_bases_aux_curves[5].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    // Compute pairs of point subtractions for later computing the pairings between them
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
        &pub_params.five_torsion_order.to_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );
    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(&prv_key.q.to_le_bytes(), prv_key.q.nbits());
    let eUV_aux_is_eUV_power_q = eUV_aux.equals(&eUV_power_q);

    let eXV = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &XV_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );
    let eXmU = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &XmU_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );

    let eYV = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &YV_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_le_bytes(),
        pub_params.five_torsion_order.nbits(),
    );
    let eYmU = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &YmU_aux_curve.to_pointx().x(),
        &pub_params.five_torsion_order.to_le_bytes(),
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
        .mul_into(&mut X_aux_curve, &U, &x.to_le_bytes(), x.nbits());
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_aux_curve, &V, &y.to_le_bytes(), y.nbits());
    ciphertext
        .shared_end_curve
        .add_into(&mut XY_aux_curve, &X_aux_curve, &Y_aux_curve);

    // [q] * X_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut U_aux_curve,
        &XY_aux_curve,
        &prv_key.q.to_le_bytes(),
        prv_key.q.nbits(),
    );
    // [2^(a-2) - q] * X_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut V_aux_curve,
        &XY_aux_curve,
        &dual_factor.to_le_bytes(),
        dual_factor.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    UV_aux_curve = U_aux_curve;
    UV_aux_curve.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

    let X_AB = ciphertext.shared_end_curve.mul(
        &UV_aux_curve,
        &prv_key.delta.to_le_bytes(),
        prv_key.delta.nbits(),
    );

    ciphertext
        .shared_end_curve
        .mul_into(&mut X_aux_curve, &U, &w.to_le_bytes(), w.nbits());
    ciphertext
        .shared_end_curve
        .mul_into(&mut Y_aux_curve, &V, &z.to_le_bytes(), z.nbits());
    ciphertext
        .shared_end_curve
        .add_into(&mut XY_aux_curve, &X_aux_curve, &Y_aux_curve);

    // [q] * Y_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut U_aux_curve,
        &XY_aux_curve,
        &prv_key.q.to_le_bytes(),
        prv_key.q.nbits(),
    );
    // [2^(a-2) - q] * Y_AB'
    ciphertext.shared_end_curve.mul_into(
        &mut V_aux_curve,
        &XY_aux_curve,
        &dual_factor.to_le_bytes(),
        dual_factor.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    UV_aux_curve = U_aux_curve;
    UV_aux_curve.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

    let Y_AB = ciphertext.shared_end_curve.mul(
        &UV_aux_curve,
        &prv_key.delta.to_le_bytes(),
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
