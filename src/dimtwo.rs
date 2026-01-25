use core::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};

use crate::{
    SUCCESS_RETVAL,
    bn::BigNum,
    dlp::{
        solve_dlp_order_powers_of_two_three, solve_dlp_order_powers_of_two_three_five,
        solve_dlp_small_prime_power_order,
    },
    endomorphism::{
        apply_endomorphism_from_quaternion, find_kernel_of_backtracking_isogeny_prime_power_degree,
        represent_integer,
    },
    inke,
    masking::mask_basis_by_same_scalar,
    poke,
    rand::sample_random_torsion_basis,
};

pub fn generate_2d_isogeny_inke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_233: usize,
>(
    pub_params: &inke::PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_233,
    >,
    q: &BigNum<NUM_WORDS_2>,
    q_dual: &BigNum<NUM_WORDS_2>,
) -> (
    EllipticProduct<Fp2>,
    (ProductPoint<Fp2>, ProductPoint<Fp2>),
    u32,
) {
    let mut retval = SUCCESS_RETVAL;

    // FIXME: When I manage to get around the shortcomings of const generics, I can use
    // widening_mul() here instead, and only progressively widen the BigNums as I multiply
    let norm: BigNum<NUM_WORDS_223> =
        q.widen() * q_dual.widen() * pub_params.three_torsion_order.widen();
    let (theta, ok) = represent_integer(&norm, &pub_params.field_characteristic);
    retval &= ok;

    /* Construct backtracking isogeny */

    // Find kernel of degree-3^b backtracking isogeny of theta
    let K3 = find_kernel_of_backtracking_isogeny_prime_power_degree(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &pub_params.three_torsion_basis,
        &pub_params.reduced_three_torsion_order,
    );

    // Apply composition of backtracking isogeny with endomorphism to the 2^a-torsion basis to obtain the kernel of the 2D-isogeny
    let (P, Q) = pub_params
        .starting_curve
        .lift_basis(&pub_params.two_torsion_basis);
    let theta_P = apply_endomorphism_from_quaternion(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &P,
    );
    let theta_Q = apply_endomorphism_from_quaternion(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &Q,
    );
    let mut theta_two_torsion_basis = [
        theta_P.to_pointx(),
        theta_Q.to_pointx(),
        pub_params
            .starting_curve
            .sub(&theta_P, &theta_Q)
            .to_pointx(),
    ];

    let (codomain, ok) = pub_params.starting_curve.three_isogeny_chain(
        &K3.to_pointx(),
        pub_params.three_torsion_exp,
        &mut theta_two_torsion_basis,
    );
    retval &= ok;

    let theta_backtrack_two_torsion_basis_img = BasisX::from_slice(&theta_two_torsion_basis);
    let (theta_backtrack_P, theta_backtrack_Q) =
        codomain.lift_basis(&theta_backtrack_two_torsion_basis_img);

    // Get rid of the 3^b-factor introduced by the backtracking isogeny
    let (scaled_theta_backtrack_P, scaled_theta_backtrack_Q) = mask_basis_by_same_scalar(
        &codomain,
        &(theta_backtrack_P, theta_backtrack_Q),
        &pub_params
            .three_torsion_order
            .invert_mod(&pub_params.full_two_torsion_order),
    );

    let two_torsion_basis_E0 = pub_params
        .starting_curve
        .lift_basis(&pub_params.two_torsion_basis);
    let (scaled_P, scaled_Q) =
        mask_basis_by_same_scalar(&pub_params.starting_curve, &two_torsion_basis_E0, &q);

    // Construct secret 2D-isogeny
    let domain = EllipticProduct::new(&pub_params.starting_curve, &codomain);
    let P1P2 = ProductPoint::new(&scaled_P, &scaled_theta_backtrack_P);
    let Q1Q2 = ProductPoint::new(&scaled_Q, &scaled_theta_backtrack_Q);

    (domain, (P1P2, Q1Q2), retval)
}

pub fn generate_2d_isogeny_poke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
>(
    pub_params: &poke::PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_23,
        NUM_WORDS_25,
        NUM_WORDS_35,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_2235,
        NUM_WORDS_2335,
        NUM_WORDS_2355,
    >,
    q: &BigNum<NUM_WORDS_2>,
    q_dual: &BigNum<NUM_WORDS_2>,
) -> (
    EllipticProduct<Fp2>,
    (ProductPoint<Fp2>, ProductPoint<Fp2>),
    u32,
) {
    let mut retval = SUCCESS_RETVAL;

    // FIXME: When I manage to get around the shortcomings of const generics, I can use
    // widening_mul() here instead, and only progressively widen the BigNums as I multiply
    let norm: BigNum<NUM_WORDS_2235> =
        q.widen() * q_dual.widen() * pub_params.three_torsion_order.widen();
    let (theta, ok) = represent_integer(&norm, &pub_params.field_characteristic);
    retval &= ok;

    /* Construct backtracking isogeny */

    // Find kernel of degree-3^b backtracking isogeny of theta
    let K3 = find_kernel_of_backtracking_isogeny_prime_power_degree(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &pub_params.three_torsion_basis,
        &pub_params.reduced_three_torsion_order,
    );

    // Apply composition of backtracking isogeny with endomorphism to the 2^a-torsion basis to obtain the kernel of the 2D-isogeny
    let (P, Q) = pub_params
        .starting_curve
        .lift_basis(&pub_params.two_torsion_basis);
    let theta_P = apply_endomorphism_from_quaternion(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &P,
    );
    let theta_Q = apply_endomorphism_from_quaternion(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &Q,
    );
    let mut theta_two_torsion_basis = [
        theta_P.to_pointx(),
        theta_Q.to_pointx(),
        pub_params
            .starting_curve
            .sub(&theta_P, &theta_Q)
            .to_pointx(),
    ];

    let (codomain, ok) = pub_params.starting_curve.three_isogeny_chain(
        &K3.to_pointx(),
        pub_params.three_torsion_exp,
        &mut theta_two_torsion_basis,
    );
    retval &= ok;

    let theta_backtrack_two_torsion_basis_img = BasisX::from_slice(&theta_two_torsion_basis);
    let (theta_backtrack_P, theta_backtrack_Q) =
        codomain.lift_basis(&theta_backtrack_two_torsion_basis_img);

    // Get rid of the 3^b-factor introduced by the backtracking isogeny
    let (scaled_theta_backtrack_P, scaled_theta_backtrack_Q) = mask_basis_by_same_scalar(
        &codomain,
        &(theta_backtrack_P, theta_backtrack_Q),
        &pub_params
            .three_torsion_order
            .invert_mod(&pub_params.full_two_torsion_order),
    );

    let two_torsion_basis_E0 = pub_params
        .starting_curve
        .lift_basis(&pub_params.two_torsion_basis);
    let (scaled_P, scaled_Q) =
        mask_basis_by_same_scalar(&pub_params.starting_curve, &two_torsion_basis_E0, &q);

    // Construct secret 2D-isogeny
    let domain = EllipticProduct::new(&pub_params.starting_curve, &codomain);
    let P1P2 = ProductPoint::new(&scaled_P, &scaled_theta_backtrack_P);
    let Q1Q2 = ProductPoint::new(&scaled_Q, &scaled_theta_backtrack_Q);

    (domain, (P1P2, Q1Q2), retval)
}

pub fn eval_2d_two_isogeny_chain_on_prime_power_torsion_basis<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_ORD: usize,
    const NUM_WORDS_COF: usize,
>(
    domain: &EllipticProduct<Fp2>,
    kernel: (&ProductPoint<Fp2>, &ProductPoint<Fp2>),
    chain_length: usize,
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
    torsion_basis: &BasisX<Fp2>,
    torsion_basis_order_base: u8,
    torsion_basis_order_exp: usize,
    torsion_basis_order: &BigNum<NUM_WORDS_ORD>,
    torsion_basis_cofactor: &BigNum<NUM_WORDS_COF>,
    p_adic_basis_for_torsion_basis_order_base: &[BigNum<NUM_WORDS_ORD>],
) -> (Curve<Fp2>, (Point<Fp2>, Point<Fp2>), u32) {
    let mut retval = SUCCESS_RETVAL;

    let (embedded_isogeny_domain, embedded_isogeny_codomain) = domain.curves();
    let (P1P2, Q1Q2) = kernel;

    // Lift basis to full points, as the 2D-isogeny function requires this as input
    let (P, Q) = embedded_isogeny_domain.lift_basis(torsion_basis);
    let PQ = embedded_isogeny_domain.sub(&P, &Q);

    // Generate random basis of the 5^c-torsion on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[torsion_basis_order_base],
        torsion_basis_order,
        torsion_basis_cofactor,
    );
    let UV = embedded_isogeny_codomain.sub(&U, &V);

    // Compute Phi' on E_domain[N] basis and random E_codomain[N] basis
    let (aux_curves, torsion_bases_on_aux_curves, ok) = domain.elliptic_product_isogeny(
        P1P2,
        Q1Q2,
        chain_length,
        &[
            ProductPoint::new(&P, &Point::INFINITY),
            ProductPoint::new(&Q, &Point::INFINITY),
            ProductPoint::new(&PQ, &Point::INFINITY),
            ProductPoint::new(&Point::INFINITY, &U),
            ProductPoint::new(&Point::INFINITY, &V),
            ProductPoint::new(&Point::INFINITY, &UV),
        ],
    );
    retval &= ok;

    let aux_curve = aux_curves.curves().0;

    // Correct the pairs of image points to overall sign
    let mut P_aux_curve = torsion_bases_on_aux_curves[0].points().0;
    let mut Q_aux_curve = torsion_bases_on_aux_curves[1].points().0;
    let mut PQ_aux_curve = torsion_bases_on_aux_curves[2].points().0;
    Q_aux_curve.set_condneg(
        !aux_curve
            .sub(&P_aux_curve, &Q_aux_curve)
            .to_pointx()
            .equals(&PQ_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = torsion_bases_on_aux_curves[3].points().0;
    let mut V_aux_curve = torsion_bases_on_aux_curves[4].points().0;
    let UV_aux_curve = torsion_bases_on_aux_curves[5].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    /* Find change-of-basis matrix */

    // Compute pairs of point subtractions for later computing the pairings between them
    let PV_aux_curve = aux_curve.sub(&P_aux_curve, &V_aux_curve);
    let PmU_aux_curve = aux_curve.add(&P_aux_curve, &U_aux_curve);

    let QV_aux_curve = aux_curve.sub(&Q_aux_curve, &V_aux_curve);
    let QmU_aux_curve = aux_curve.add(&Q_aux_curve, &U_aux_curve);

    // Compute the pairings e(U, V), e(X, V) = e(U, V)^x and e(X, -U) = e(U, V)^y,
    // e(Y, V) = e(U, V)^w and e(Y, -U) = e(U, V)^z
    let eUV_aux = aux_curve.weil_pairing(
        &U_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &UV_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let ePV = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &PV_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let ePmU = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &PmU_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let eQV = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &QV_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let eQmU = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &QmU_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );

    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(&degree.to_le_bytes(), degree.nbits());
    let dual_has_degree_q = eUV_aux.equals(&eUV_power_q);

    // Solve discrete logarithm between pairings to obtain expression of P' in terms of <U',V'>
    let (x, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &ePV,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;
    let (y, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &ePmU,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;
    let (w, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eQV,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;
    let (z, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eQmU,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;

    /* Derive final image points by applying the change-of-basis matrix */

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &x.to_le_bytes(), x.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &y.to_le_bytes(), y.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_P = U_aux_curve;
    img_P.set_cond(&V_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &w.to_le_bytes(), w.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &z.to_le_bytes(), z.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Q = U_aux_curve;
    img_Q.set_cond(&V_aux_curve, !dual_has_degree_q);

    (aux_curve, (img_P, img_Q), retval)
}

pub fn eval_2d_two_isogeny_chain_inke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_233: usize,
>(
    domain: &EllipticProduct<Fp2>,
    kernel: (&ProductPoint<Fp2>, &ProductPoint<Fp2>),
    chain_length: usize,
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
    full_torsion_basis: &BasisX<Fp2>,
    torsion_basis_prime_power_exps: (usize, usize),
    torsion_basis_orders: (&BigNum<NUM_WORDS_2>, &BigNum<NUM_WORDS_3>),
    full_torsion_basis_order: &BigNum<NUM_WORDS_P>,
    full_torsion_basis_cofactor: &BigNum<1>,
    p_adic_bases: (&[BigNum<NUM_WORDS_2>], &[BigNum<NUM_WORDS_3>]),
    intermediate_bignum_sizes: PhantomData<([(); NUM_WORDS_223], [(); NUM_WORDS_233])>,
) -> (
    (Point<Fp2>, Point<Fp2>),
    Curve<Fp2>,
    (Point<Fp2>, Point<Fp2>),
    u32,
) {
    let mut retval = SUCCESS_RETVAL;

    let (embedded_isogeny_domain, embedded_isogeny_codomain) = domain.curves();
    let (P1P2, Q1Q2) = kernel;

    // Lift basis to full points, as the 2D-isogeny function requires this as input
    let (P, Q) = embedded_isogeny_domain.lift_basis(full_torsion_basis);
    let PQ = embedded_isogeny_domain.sub(&P, &Q);

    // Generate random basis of the (2^a * 3^b)-torsion on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[2, 3],
        full_torsion_basis_order,
        full_torsion_basis_cofactor,
    );
    let UV = embedded_isogeny_codomain.sub(&U, &V);

    // Compute Phi' on E_domain[N] basis and random E_codomain[N] basis
    let (aux_curves, torsion_bases_on_aux_curves, ok) = domain.elliptic_product_isogeny(
        P1P2,
        Q1Q2,
        chain_length,
        &[
            ProductPoint::new(&P, &Point::INFINITY),
            ProductPoint::new(&Q, &Point::INFINITY),
            ProductPoint::new(&PQ, &Point::INFINITY),
            ProductPoint::new(&Point::INFINITY, &U),
            ProductPoint::new(&Point::INFINITY, &V),
            ProductPoint::new(&Point::INFINITY, &UV),
        ],
    );
    retval &= ok;

    let aux_curve = aux_curves.curves().0;

    // Correct the pairs of image points to overall sign
    let mut P_aux_curve = torsion_bases_on_aux_curves[0].points().0;
    let mut Q_aux_curve = torsion_bases_on_aux_curves[1].points().0;
    let mut PQ_aux_curve = torsion_bases_on_aux_curves[2].points().0;
    Q_aux_curve.set_condneg(
        !aux_curve
            .sub(&P_aux_curve, &Q_aux_curve)
            .to_pointx()
            .equals(&PQ_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = torsion_bases_on_aux_curves[3].points().0;
    let mut V_aux_curve = torsion_bases_on_aux_curves[4].points().0;
    let UV_aux_curve = torsion_bases_on_aux_curves[5].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    /* Find change-of-basis matrix */

    // Compute pairs of point subtractions for later computing the pairings between them
    let PV_aux_curve = aux_curve.sub(&P_aux_curve, &V_aux_curve);
    let PmU_aux_curve = aux_curve.add(&P_aux_curve, &U_aux_curve);

    let QV_aux_curve = aux_curve.sub(&Q_aux_curve, &V_aux_curve);
    let QmU_aux_curve = aux_curve.add(&Q_aux_curve, &U_aux_curve);

    // Compute the pairings e(U, V), e(X, V) = e(U, V)^x and e(X, -U) = e(U, V)^y,
    // e(Y, V) = e(U, V)^w and e(Y, -U) = e(U, V)^z
    let eUV_aux = aux_curve.weil_pairing(
        &U_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &UV_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let ePV = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &PV_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let ePmU = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &PmU_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let eQV = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &QV_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let eQmU = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &QmU_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );

    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(&degree.to_le_bytes(), degree.nbits());
    let dual_has_degree_q = eUV_aux.equals(&eUV_power_q);

    // Solve discrete logarithm between pairings to obtain expression of P' in terms of <U',V'>
    let (x, ok) = solve_dlp_order_powers_of_two_three(
        &eUV_aux,
        &ePV,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;
    let (y, ok) = solve_dlp_order_powers_of_two_three(
        &eUV_aux,
        &ePmU,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;
    let (w, ok) = solve_dlp_order_powers_of_two_three(
        &eUV_aux,
        &eQV,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;
    let (z, ok) = solve_dlp_order_powers_of_two_three(
        &eUV_aux,
        &eQmU,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;

    /* Derive final image points on intermediate curve */

    let img_P_intermediate = P_aux_curve.clone();
    let img_Q_intermediate = Q_aux_curve.clone();

    /* Derive final image points on embedded codomain curve by applying the change-of-basis matrix */

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &x.to_le_bytes(), x.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &y.to_le_bytes(), y.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_P_codomain = U_aux_curve;
    img_P_codomain.set_cond(&V_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &w.to_le_bytes(), w.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &z.to_le_bytes(), z.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Q_codomain = U_aux_curve;
    img_Q_codomain.set_cond(&V_aux_curve, !dual_has_degree_q);

    (
        (img_P_codomain, img_Q_codomain),
        aux_curve,
        (img_P_intermediate, img_Q_intermediate),
        retval,
    )
}

pub fn eval_2d_two_isogeny_chain_inke_separate_bases<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
>(
    domain: &EllipticProduct<Fp2>,
    kernel: (&ProductPoint<Fp2>, &ProductPoint<Fp2>),
    chain_length: usize,
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
    torsion_bases: (&BasisX<Fp2>, &BasisX<Fp2>),
    torsion_basis_prime_power_exps: (usize, usize),
    torsion_basis_orders: (&BigNum<NUM_WORDS_2>, &BigNum<NUM_WORDS_3>),
    torsion_basis_cofactors: (&BigNum<NUM_WORDS_3>, &BigNum<NUM_WORDS_2>),
    p_adic_bases: (&[BigNum<NUM_WORDS_2>], &[BigNum<NUM_WORDS_3>]),
) -> (
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    Curve<Fp2>,
    (Point<Fp2>, Point<Fp2>),
    u32,
) {
    let mut retval = SUCCESS_RETVAL;

    let (embedded_isogeny_domain, embedded_isogeny_codomain) = domain.curves();
    let (P1P2, Q1Q2) = kernel;

    // Lift bases to full points, as the 2D-isogeny function requires these as input
    let (P, Q) = embedded_isogeny_domain.lift_basis(torsion_bases.0);
    let PQ = embedded_isogeny_domain.sub(&P, &Q);
    let (R, S) = embedded_isogeny_domain.lift_basis(torsion_bases.1);
    let RS = embedded_isogeny_domain.sub(&R, &S);

    // Generate random bases of the different torsion subgroups on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[2],
        torsion_basis_orders.0,
        torsion_basis_cofactors.0,
    );
    let UV = embedded_isogeny_codomain.sub(&U, &V);
    let (W, Z, _) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[3],
        torsion_basis_orders.1,
        torsion_basis_cofactors.1,
    );
    let WZ = embedded_isogeny_codomain.sub(&W, &Z);

    // Compute Phi' on E_domain[N] basis and random E_codomain[N] basis
    let (aux_curves, torsion_bases_on_aux_curves, ok) = domain.elliptic_product_isogeny(
        P1P2,
        Q1Q2,
        chain_length,
        &[
            ProductPoint::new(&P, &Point::INFINITY),
            ProductPoint::new(&Q, &Point::INFINITY),
            ProductPoint::new(&PQ, &Point::INFINITY),
            ProductPoint::new(&R, &Point::INFINITY),
            ProductPoint::new(&S, &Point::INFINITY),
            ProductPoint::new(&RS, &Point::INFINITY),
            ProductPoint::new(&Point::INFINITY, &U),
            ProductPoint::new(&Point::INFINITY, &V),
            ProductPoint::new(&Point::INFINITY, &UV),
            ProductPoint::new(&Point::INFINITY, &W),
            ProductPoint::new(&Point::INFINITY, &Z),
            ProductPoint::new(&Point::INFINITY, &WZ),
        ],
    );
    retval &= ok;

    let aux_curve = aux_curves.curves().0;

    // Correct the pairs of image points to overall sign
    let mut P_aux_curve = torsion_bases_on_aux_curves[0].points().0;
    let mut Q_aux_curve = torsion_bases_on_aux_curves[1].points().0;
    let mut PQ_aux_curve = torsion_bases_on_aux_curves[2].points().0;
    Q_aux_curve.set_condneg(
        !aux_curve
            .sub(&P_aux_curve, &Q_aux_curve)
            .to_pointx()
            .equals(&PQ_aux_curve.to_pointx()),
    );

    let mut R_aux_curve = torsion_bases_on_aux_curves[3].points().0;
    let mut S_aux_curve = torsion_bases_on_aux_curves[4].points().0;
    let mut RS_aux_curve = torsion_bases_on_aux_curves[5].points().0;
    S_aux_curve.set_condneg(
        !aux_curve
            .sub(&R_aux_curve, &S_aux_curve)
            .to_pointx()
            .equals(&RS_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = torsion_bases_on_aux_curves[6].points().0;
    let mut V_aux_curve = torsion_bases_on_aux_curves[7].points().0;
    let UV_aux_curve = torsion_bases_on_aux_curves[8].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    let mut W_aux_curve = torsion_bases_on_aux_curves[9].points().0;
    let mut Z_aux_curve = torsion_bases_on_aux_curves[10].points().0;
    let WZ_aux_curve = torsion_bases_on_aux_curves[11].points().0;
    Z_aux_curve.set_condneg(
        !aux_curve
            .sub(&W_aux_curve, &Z_aux_curve)
            .to_pointx()
            .equals(&WZ_aux_curve.to_pointx()),
    );

    /* Find change-of-basis matrix */

    // Compute pairs of point subtractions for later computing the pairings between them
    let PV_aux_curve = aux_curve.sub(&P_aux_curve, &V_aux_curve);
    let PmU_aux_curve = aux_curve.add(&P_aux_curve, &U_aux_curve);
    let QV_aux_curve = aux_curve.sub(&Q_aux_curve, &V_aux_curve);
    let QmU_aux_curve = aux_curve.add(&Q_aux_curve, &U_aux_curve);

    let RZ_aux_curve = aux_curve.sub(&R_aux_curve, &Z_aux_curve);
    let RmW_aux_curve = aux_curve.add(&R_aux_curve, &W_aux_curve);
    let SZ_aux_curve = aux_curve.sub(&S_aux_curve, &Z_aux_curve);
    let SmW_aux_curve = aux_curve.add(&S_aux_curve, &W_aux_curve);

    // Compute the pairings between respective pairs of input and random torsion bases,
    // which will be the generators of the discrete log subgroups we attempt to solve
    let eUV_aux = aux_curve.weil_pairing(
        &U_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &UV_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let eWZ_aux = aux_curve.weil_pairing(
        &W_aux_curve.to_pointx().x(),
        &Z_aux_curve.to_pointx().x(),
        &WZ_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );

    // Compute the pairings for which we will solve the discrete log w.r.t.
    // the above generators for the entries of the change-of-basis matrices
    let ePV = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &PV_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let ePmU = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &PmU_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let eQV = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &QV_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let eQmU = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &QmU_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );

    let eRZ = aux_curve.weil_pairing(
        &R_aux_curve.to_pointx().x(),
        &Z_aux_curve.to_pointx().x(),
        &RZ_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eRmW = aux_curve.weil_pairing(
        &R_aux_curve.to_pointx().x(),
        &W_aux_curve.to_pointx().x(),
        &RmW_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eSZ = aux_curve.weil_pairing(
        &S_aux_curve.to_pointx().x(),
        &Z_aux_curve.to_pointx().x(),
        &SZ_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eSmW = aux_curve.weil_pairing(
        &S_aux_curve.to_pointx().x(),
        &W_aux_curve.to_pointx().x(),
        &SmW_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );

    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(&degree.to_le_bytes(), degree.nbits());
    let dual_has_degree_q = eUV_aux.equals(&eUV_power_q);

    // Solve discrete logarithm between pairings to obtain expression of P',Q' in terms of <U',V'>
    let (x2, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &ePV,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;
    let (y2, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &ePmU,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;
    let (w2, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eQV,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;
    let (z2, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eQmU,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;

    // Solve discrete logarithm between pairings to obtain expression of R',S' in terms of <W',Z'>
    let (x3, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eRZ,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;
    let (y3, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eRmW,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;
    let (w3, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eSZ,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;
    let (z3, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eSmW,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;

    /* Derive final image points on intermediate curve */

    let img_R_intermediate = R_aux_curve.clone();
    let img_S_intermediate = S_aux_curve.clone();

    /* Derive final image points on embedded codomain curve by applying the change-of-basis matrix */

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &x2.to_le_bytes(), x2.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &y2.to_le_bytes(), y2.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_P_codomain = U_aux_curve;
    img_P_codomain.set_cond(&V_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &w2.to_le_bytes(), w2.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &z2.to_le_bytes(), z2.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Q_codomain = U_aux_curve;
    img_Q_codomain.set_cond(&V_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut R_aux_curve, &W, &x3.to_le_bytes(), x3.nbits());
    embedded_isogeny_codomain.mul_into(&mut S_aux_curve, &Z, &y3.to_le_bytes(), y3.nbits());
    embedded_isogeny_codomain.add_into(&mut RS_aux_curve, &R_aux_curve, &S_aux_curve);

    // [q] * R'
    embedded_isogeny_codomain.mul_into(
        &mut W_aux_curve,
        &RS_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * S'
    embedded_isogeny_codomain.mul_into(
        &mut Z_aux_curve,
        &RS_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_R_codomain = W_aux_curve;
    img_R_codomain.set_cond(&Z_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut R_aux_curve, &W, &w3.to_le_bytes(), w3.nbits());
    embedded_isogeny_codomain.mul_into(&mut S_aux_curve, &Z, &z3.to_le_bytes(), z3.nbits());
    embedded_isogeny_codomain.add_into(&mut RS_aux_curve, &R_aux_curve, &S_aux_curve);

    // [q] * S'
    embedded_isogeny_codomain.mul_into(
        &mut W_aux_curve,
        &RS_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * S'
    embedded_isogeny_codomain.mul_into(
        &mut Z_aux_curve,
        &RS_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_S_codomain = W_aux_curve;
    img_S_codomain.set_cond(&Z_aux_curve, !dual_has_degree_q);

    (
        (img_P_codomain, img_Q_codomain),
        (img_R_codomain, img_S_codomain),
        aux_curve,
        (img_R_intermediate, img_S_intermediate),
        retval,
    )
}

pub fn eval_2d_two_isogeny_chain_poke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_235: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
>(
    domain: &EllipticProduct<Fp2>,
    kernel: (&ProductPoint<Fp2>, &ProductPoint<Fp2>),
    chain_length: usize,
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
    full_torsion_basis: &BasisX<Fp2>,
    torsion_basis_prime_power_exps: (usize, usize, usize),
    torsion_basis_orders: (
        &BigNum<NUM_WORDS_2>,
        &BigNum<NUM_WORDS_3>,
        &BigNum<NUM_WORDS_5>,
    ),
    partial_products_of_prime_powers: (
        &BigNum<NUM_WORDS_35>,
        &BigNum<NUM_WORDS_25>,
        &BigNum<NUM_WORDS_23>,
    ),
    full_torsion_basis_order: &BigNum<NUM_WORDS_235>,
    full_torsion_basis_cofactor: &BigNum<1>,
    p_adic_bases: (
        &[BigNum<NUM_WORDS_2>],
        &[BigNum<NUM_WORDS_3>],
        &[BigNum<NUM_WORDS_5>],
    ),
    intermediate_bignum_sizes: PhantomData<(
        [(); NUM_WORDS_2235],
        [(); NUM_WORDS_2335],
        [(); NUM_WORDS_2355],
    )>,
) -> ((Point<Fp2>, Point<Fp2>), u32) {
    let mut retval = SUCCESS_RETVAL;

    let (embedded_isogeny_domain, embedded_isogeny_codomain) = domain.curves();
    let (P1P2, Q1Q2) = kernel;

    // Lift basis to full points, as the 2D-isogeny function requires this as input
    let (P, Q) = embedded_isogeny_domain.lift_basis(full_torsion_basis);
    let PQ = embedded_isogeny_domain.sub(&P, &Q);

    // Generate random basis of the (2^a * 3^b * 5^c)-torsion on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[2, 3, 5],
        full_torsion_basis_order,
        full_torsion_basis_cofactor,
    );
    let UV = embedded_isogeny_codomain.sub(&U, &V);

    // Compute Phi' on E_domain[N] basis and random E_codomain[N] basis
    let (aux_curves, torsion_bases_on_aux_curves, ok) = domain.elliptic_product_isogeny(
        P1P2,
        Q1Q2,
        chain_length,
        &[
            ProductPoint::new(&P, &Point::INFINITY),
            ProductPoint::new(&Q, &Point::INFINITY),
            ProductPoint::new(&PQ, &Point::INFINITY),
            ProductPoint::new(&Point::INFINITY, &U),
            ProductPoint::new(&Point::INFINITY, &V),
            ProductPoint::new(&Point::INFINITY, &UV),
        ],
    );
    retval &= ok;

    let aux_curve = aux_curves.curves().0;

    // Correct the pairs of image points to overall sign
    let mut P_aux_curve = torsion_bases_on_aux_curves[0].points().0;
    let mut Q_aux_curve = torsion_bases_on_aux_curves[1].points().0;
    let mut PQ_aux_curve = torsion_bases_on_aux_curves[2].points().0;
    Q_aux_curve.set_condneg(
        !aux_curve
            .sub(&P_aux_curve, &Q_aux_curve)
            .to_pointx()
            .equals(&PQ_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = torsion_bases_on_aux_curves[3].points().0;
    let mut V_aux_curve = torsion_bases_on_aux_curves[4].points().0;
    let UV_aux_curve = torsion_bases_on_aux_curves[5].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    /* Find change-of-basis matrix */

    // Compute pairs of point subtractions for later computing the pairings between them
    let PV_aux_curve = aux_curve.sub(&P_aux_curve, &V_aux_curve);
    let PmU_aux_curve = aux_curve.add(&P_aux_curve, &U_aux_curve);

    let QV_aux_curve = aux_curve.sub(&Q_aux_curve, &V_aux_curve);
    let QmU_aux_curve = aux_curve.add(&Q_aux_curve, &U_aux_curve);

    // Compute the pairings e(U, V), e(X, V) = e(U, V)^x and e(X, -U) = e(U, V)^y,
    // e(Y, V) = e(U, V)^w and e(Y, -U) = e(U, V)^z
    let eUV_aux = aux_curve.weil_pairing(
        &U_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &UV_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let ePV = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &PV_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let ePmU = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &PmU_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let eQV = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &QV_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );
    let eQmU = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &QmU_aux_curve.to_pointx().x(),
        &full_torsion_basis_order.to_le_bytes(),
        full_torsion_basis_order.nbits(),
    );

    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(&degree.to_le_bytes(), degree.nbits());
    let dual_has_degree_q = eUV_aux.equals(&eUV_power_q);

    // Solve discrete logarithm between pairings to obtain expression of P' in terms of <U',V'>
    let (x, ok) = solve_dlp_order_powers_of_two_three_five(
        &eUV_aux,
        &ePV,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        partial_products_of_prime_powers,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;
    let (y, ok) = solve_dlp_order_powers_of_two_three_five(
        &eUV_aux,
        &ePmU,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        partial_products_of_prime_powers,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;
    let (w, ok) = solve_dlp_order_powers_of_two_three_five(
        &eUV_aux,
        &eQV,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        partial_products_of_prime_powers,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;
    let (z, ok) = solve_dlp_order_powers_of_two_three_five(
        &eUV_aux,
        &eQmU,
        torsion_basis_prime_power_exps,
        torsion_basis_orders,
        partial_products_of_prime_powers,
        full_torsion_basis_order,
        p_adic_bases,
        intermediate_bignum_sizes,
    );
    retval &= ok;

    /* Derive final image points by applying the change-of-basis matrix */

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &x.to_le_bytes(), x.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &y.to_le_bytes(), y.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_P = U_aux_curve;
    img_P.set_cond(&V_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrix to points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &U, &w.to_le_bytes(), w.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &V, &z.to_le_bytes(), z.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Q = U_aux_curve;
    img_Q.set_cond(&V_aux_curve, !dual_has_degree_q);

    ((img_P, img_Q), retval)
}

pub fn eval_2d_two_isogeny_chain_poke_separate_bases<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
>(
    domain: &EllipticProduct<Fp2>,
    kernel: (&ProductPoint<Fp2>, &ProductPoint<Fp2>),
    chain_length: usize,
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
    torsion_bases: (&BasisX<Fp2>, &BasisX<Fp2>, &BasisX<Fp2>),
    torsion_basis_prime_power_exps: (usize, usize, usize),
    torsion_basis_orders: (
        &BigNum<NUM_WORDS_2>,
        &BigNum<NUM_WORDS_3>,
        &BigNum<NUM_WORDS_5>,
    ),
    torsion_basis_cofactors: (
        &BigNum<NUM_WORDS_35>,
        &BigNum<NUM_WORDS_25>,
        &BigNum<NUM_WORDS_23>,
    ),
    p_adic_bases: (
        &[BigNum<NUM_WORDS_2>],
        &[BigNum<NUM_WORDS_3>],
        &[BigNum<NUM_WORDS_5>],
    ),
) -> (
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    u32,
) {
    let mut retval = SUCCESS_RETVAL;

    let (embedded_isogeny_domain, embedded_isogeny_codomain) = domain.curves();
    let (P1P2, Q1Q2) = kernel;

    // Lift bases to full points, as the 2D-isogeny function requires this as input
    let (P, Q) = embedded_isogeny_domain.lift_basis(torsion_bases.0);
    let PQ = embedded_isogeny_domain.sub(&P, &Q);
    let (R, S) = embedded_isogeny_domain.lift_basis(torsion_bases.1);
    let RS = embedded_isogeny_domain.sub(&R, &S);
    let (X, Y) = embedded_isogeny_domain.lift_basis(torsion_bases.2);
    let XY = embedded_isogeny_domain.sub(&X, &Y);

    // Generate random bases of the different torsion subgroups on E_AB
    let (M, N, eMN_AB) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[2],
        torsion_basis_orders.0,
        torsion_basis_cofactors.0,
    );
    let MN = embedded_isogeny_codomain.sub(&M, &N);
    let (U, V, _) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[3],
        torsion_basis_orders.1,
        torsion_basis_cofactors.1,
    );
    let UV = embedded_isogeny_codomain.sub(&U, &V);
    let (W, Z, _) = sample_random_torsion_basis(
        &embedded_isogeny_codomain,
        &[5],
        torsion_basis_orders.2,
        torsion_basis_cofactors.2,
    );
    let WZ = embedded_isogeny_codomain.sub(&W, &Z);

    // Compute Phi' on E_domain[N] basis and random E_codomain[N] basis
    let (aux_curves, torsion_bases_on_aux_curves, ok) = domain.elliptic_product_isogeny(
        P1P2,
        Q1Q2,
        chain_length,
        &[
            ProductPoint::new(&P, &Point::INFINITY),
            ProductPoint::new(&Q, &Point::INFINITY),
            ProductPoint::new(&PQ, &Point::INFINITY),
            ProductPoint::new(&R, &Point::INFINITY),
            ProductPoint::new(&S, &Point::INFINITY),
            ProductPoint::new(&RS, &Point::INFINITY),
            ProductPoint::new(&X, &Point::INFINITY),
            ProductPoint::new(&Y, &Point::INFINITY),
            ProductPoint::new(&XY, &Point::INFINITY),
            ProductPoint::new(&Point::INFINITY, &M),
            ProductPoint::new(&Point::INFINITY, &N),
            ProductPoint::new(&Point::INFINITY, &MN),
            ProductPoint::new(&Point::INFINITY, &U),
            ProductPoint::new(&Point::INFINITY, &V),
            ProductPoint::new(&Point::INFINITY, &UV),
            ProductPoint::new(&Point::INFINITY, &W),
            ProductPoint::new(&Point::INFINITY, &Z),
            ProductPoint::new(&Point::INFINITY, &WZ),
        ],
    );
    retval &= ok;

    let aux_curve = aux_curves.curves().0;

    // Correct the pairs of image points to overall sign
    let mut P_aux_curve = torsion_bases_on_aux_curves[0].points().0;
    let mut Q_aux_curve = torsion_bases_on_aux_curves[1].points().0;
    let mut PQ_aux_curve = torsion_bases_on_aux_curves[2].points().0;
    Q_aux_curve.set_condneg(
        !aux_curve
            .sub(&P_aux_curve, &Q_aux_curve)
            .to_pointx()
            .equals(&PQ_aux_curve.to_pointx()),
    );

    let mut R_aux_curve = torsion_bases_on_aux_curves[3].points().0;
    let mut S_aux_curve = torsion_bases_on_aux_curves[4].points().0;
    let mut RS_aux_curve = torsion_bases_on_aux_curves[5].points().0;
    S_aux_curve.set_condneg(
        !aux_curve
            .sub(&R_aux_curve, &S_aux_curve)
            .to_pointx()
            .equals(&RS_aux_curve.to_pointx()),
    );

    let mut X_aux_curve = torsion_bases_on_aux_curves[6].points().0;
    let mut Y_aux_curve = torsion_bases_on_aux_curves[7].points().0;
    let mut XY_aux_curve = torsion_bases_on_aux_curves[8].points().0;
    Y_aux_curve.set_condneg(
        !aux_curve
            .sub(&X_aux_curve, &Y_aux_curve)
            .to_pointx()
            .equals(&XY_aux_curve.to_pointx()),
    );

    let mut M_aux_curve = torsion_bases_on_aux_curves[9].points().0;
    let mut N_aux_curve = torsion_bases_on_aux_curves[10].points().0;
    let MN_aux_curve = torsion_bases_on_aux_curves[11].points().0;
    N_aux_curve.set_condneg(
        !aux_curve
            .sub(&M_aux_curve, &N_aux_curve)
            .to_pointx()
            .equals(&MN_aux_curve.to_pointx()),
    );

    let mut U_aux_curve = torsion_bases_on_aux_curves[12].points().0;
    let mut V_aux_curve = torsion_bases_on_aux_curves[13].points().0;
    let UV_aux_curve = torsion_bases_on_aux_curves[14].points().0;
    V_aux_curve.set_condneg(
        !aux_curve
            .sub(&U_aux_curve, &V_aux_curve)
            .to_pointx()
            .equals(&UV_aux_curve.to_pointx()),
    );

    let mut W_aux_curve = torsion_bases_on_aux_curves[15].points().0;
    let mut Z_aux_curve = torsion_bases_on_aux_curves[16].points().0;
    let WZ_aux_curve = torsion_bases_on_aux_curves[17].points().0;
    Z_aux_curve.set_condneg(
        !aux_curve
            .sub(&W_aux_curve, &Z_aux_curve)
            .to_pointx()
            .equals(&WZ_aux_curve.to_pointx()),
    );

    /* Find change-of-basis matrices for each pair of basis points and its random torsion basis */

    // Compute pairs of point subtractions for later computing the pairings between them
    let PN_aux_curve = aux_curve.sub(&P_aux_curve, &N_aux_curve);
    let PmM_aux_curve = aux_curve.add(&P_aux_curve, &M_aux_curve);
    let QN_aux_curve = aux_curve.sub(&Q_aux_curve, &N_aux_curve);
    let QmM_aux_curve = aux_curve.add(&Q_aux_curve, &M_aux_curve);

    let RV_aux_curve = aux_curve.sub(&R_aux_curve, &V_aux_curve);
    let RmU_aux_curve = aux_curve.add(&R_aux_curve, &U_aux_curve);
    let SV_aux_curve = aux_curve.sub(&S_aux_curve, &V_aux_curve);
    let SmU_aux_curve = aux_curve.add(&S_aux_curve, &U_aux_curve);

    let XZ_aux_curve = aux_curve.sub(&X_aux_curve, &Z_aux_curve);
    let XmW_aux_curve = aux_curve.add(&X_aux_curve, &W_aux_curve);
    let YZ_aux_curve = aux_curve.sub(&Y_aux_curve, &Z_aux_curve);
    let YmW_aux_curve = aux_curve.add(&Y_aux_curve, &W_aux_curve);

    // Compute the pairings between respective pairs of input and random torsion bases,
    // which will be the generators of the discrete log subgroups we attempt to solve
    let eMN_aux = aux_curve.weil_pairing(
        &M_aux_curve.to_pointx().x(),
        &N_aux_curve.to_pointx().x(),
        &MN_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let eUV_aux = aux_curve.weil_pairing(
        &U_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &UV_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eWZ_aux = aux_curve.weil_pairing(
        &W_aux_curve.to_pointx().x(),
        &Z_aux_curve.to_pointx().x(),
        &WZ_aux_curve.to_pointx().x(),
        &torsion_basis_orders.2.to_le_bytes(),
        torsion_basis_orders.2.nbits(),
    );

    // Compute the pairings for which we will solve the discrete log w.r.t.
    // the above generators for the entries of the change-of-basis matrices
    let ePN = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &N_aux_curve.to_pointx().x(),
        &PN_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let ePmM = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &M_aux_curve.to_pointx().x(),
        &PmM_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let eQN = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &N_aux_curve.to_pointx().x(),
        &QN_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );
    let eQmM = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &M_aux_curve.to_pointx().x(),
        &QmM_aux_curve.to_pointx().x(),
        &torsion_basis_orders.0.to_le_bytes(),
        torsion_basis_orders.0.nbits(),
    );

    let eRV = aux_curve.weil_pairing(
        &R_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &RV_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eRmU = aux_curve.weil_pairing(
        &R_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &RmU_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eSV = aux_curve.weil_pairing(
        &S_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &SV_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );
    let eSmU = aux_curve.weil_pairing(
        &S_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &SmU_aux_curve.to_pointx().x(),
        &torsion_basis_orders.1.to_le_bytes(),
        torsion_basis_orders.1.nbits(),
    );

    let eXZ = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &Z_aux_curve.to_pointx().x(),
        &XZ_aux_curve.to_pointx().x(),
        &torsion_basis_orders.2.to_le_bytes(),
        torsion_basis_orders.2.nbits(),
    );
    let eXmW = aux_curve.weil_pairing(
        &X_aux_curve.to_pointx().x(),
        &W_aux_curve.to_pointx().x(),
        &XmW_aux_curve.to_pointx().x(),
        &torsion_basis_orders.2.to_le_bytes(),
        torsion_basis_orders.2.nbits(),
    );
    let eYZ = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &Z_aux_curve.to_pointx().x(),
        &YZ_aux_curve.to_pointx().x(),
        &torsion_basis_orders.2.to_le_bytes(),
        torsion_basis_orders.2.nbits(),
    );
    let eYmW = aux_curve.weil_pairing(
        &Y_aux_curve.to_pointx().x(),
        &W_aux_curve.to_pointx().x(),
        &YmW_aux_curve.to_pointx().x(),
        &torsion_basis_orders.2.to_le_bytes(),
        torsion_basis_orders.2.nbits(),
    );

    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eMN_power_q = eMN_AB.pow(&degree.to_le_bytes(), degree.nbits());
    let dual_has_degree_q = eMN_aux.equals(&eMN_power_q);

    // Solve discrete logarithm between pairings to obtain expression of each
    // image point in terms of the image points of its respective random torsion basis
    let (x2, ok) = solve_dlp_small_prime_power_order(
        &eMN_aux,
        &ePN,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;
    let (y2, ok) = solve_dlp_small_prime_power_order(
        &eMN_aux,
        &ePmM,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;
    let (w2, ok) = solve_dlp_small_prime_power_order(
        &eMN_aux,
        &eQN,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;
    let (z2, ok) = solve_dlp_small_prime_power_order(
        &eMN_aux,
        &eQmM,
        2,
        torsion_basis_prime_power_exps.0,
        p_adic_bases.0,
    );
    retval &= ok;

    let (x3, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eRV,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;
    let (y3, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eRmU,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;
    let (w3, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eSV,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;
    let (z3, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eSmU,
        3,
        torsion_basis_prime_power_exps.1,
        p_adic_bases.1,
    );
    retval &= ok;

    let (x5, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eXZ,
        5,
        torsion_basis_prime_power_exps.2,
        p_adic_bases.2,
    );
    retval &= ok;
    let (y5, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eXmW,
        5,
        torsion_basis_prime_power_exps.2,
        p_adic_bases.2,
    );
    retval &= ok;
    let (w5, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eYZ,
        5,
        torsion_basis_prime_power_exps.2,
        p_adic_bases.2,
    );
    retval &= ok;
    let (z5, ok) = solve_dlp_small_prime_power_order(
        &eWZ_aux,
        &eYmW,
        5,
        torsion_basis_prime_power_exps.2,
        p_adic_bases.2,
    );
    retval &= ok;

    /* Derive final image points by applying the change-of-basis matrices */

    // Apply change-of-basis matrices to respective basis points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &M, &x2.to_le_bytes(), x2.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &N, &y2.to_le_bytes(), y2.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut M_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut N_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_P = M_aux_curve;
    img_P.set_cond(&N_aux_curve, !dual_has_degree_q);

    embedded_isogeny_codomain.mul_into(&mut P_aux_curve, &M, &w2.to_le_bytes(), w2.nbits());
    embedded_isogeny_codomain.mul_into(&mut Q_aux_curve, &N, &z2.to_le_bytes(), z2.nbits());
    embedded_isogeny_codomain.add_into(&mut PQ_aux_curve, &P_aux_curve, &Q_aux_curve);

    // [q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut M_aux_curve,
        &PQ_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut N_aux_curve,
        &PQ_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Q = M_aux_curve;
    img_Q.set_cond(&N_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrices to respective basis points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut R_aux_curve, &U, &x3.to_le_bytes(), x3.nbits());
    embedded_isogeny_codomain.mul_into(&mut S_aux_curve, &V, &y3.to_le_bytes(), y3.nbits());
    embedded_isogeny_codomain.add_into(&mut RS_aux_curve, &R_aux_curve, &S_aux_curve);

    // [q] * R'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &RS_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * R'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &RS_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_R = U_aux_curve;
    img_R.set_cond(&V_aux_curve, !dual_has_degree_q);

    embedded_isogeny_codomain.mul_into(&mut R_aux_curve, &U, &w3.to_le_bytes(), w3.nbits());
    embedded_isogeny_codomain.mul_into(&mut S_aux_curve, &V, &z3.to_le_bytes(), z3.nbits());
    embedded_isogeny_codomain.add_into(&mut RS_aux_curve, &R_aux_curve, &S_aux_curve);

    // [q] * S'
    embedded_isogeny_codomain.mul_into(
        &mut U_aux_curve,
        &RS_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * S'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &RS_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_S = U_aux_curve;
    img_S.set_cond(&V_aux_curve, !dual_has_degree_q);

    // Apply change-of-basis matrices to respective basis points on embedded isogeny codomain curve
    // (reusing temporary intermediate curve points as an optimization)
    embedded_isogeny_codomain.mul_into(&mut X_aux_curve, &W, &x5.to_le_bytes(), x5.nbits());
    embedded_isogeny_codomain.mul_into(&mut Y_aux_curve, &Z, &y5.to_le_bytes(), y5.nbits());
    embedded_isogeny_codomain.add_into(&mut XY_aux_curve, &X_aux_curve, &Y_aux_curve);

    // [q] * X'
    embedded_isogeny_codomain.mul_into(
        &mut W_aux_curve,
        &XY_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * X'
    embedded_isogeny_codomain.mul_into(
        &mut Z_aux_curve,
        &XY_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_X = W_aux_curve;
    img_X.set_cond(&Z_aux_curve, !dual_has_degree_q);

    embedded_isogeny_codomain.mul_into(&mut X_aux_curve, &W, &w5.to_le_bytes(), w5.nbits());
    embedded_isogeny_codomain.mul_into(&mut Y_aux_curve, &Z, &z5.to_le_bytes(), z5.nbits());
    embedded_isogeny_codomain.add_into(&mut XY_aux_curve, &X_aux_curve, &Y_aux_curve);

    // [q] * Y'
    embedded_isogeny_codomain.mul_into(
        &mut W_aux_curve,
        &XY_aux_curve,
        &degree.to_le_bytes(),
        degree.nbits(),
    );
    // [2^(a-2) - q] * Y'
    embedded_isogeny_codomain.mul_into(
        &mut Z_aux_curve,
        &XY_aux_curve,
        &degree_dual.to_le_bytes(),
        degree_dual.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Y = W_aux_curve;
    img_Y.set_cond(&Z_aux_curve, !dual_has_degree_q);

    ((img_P, img_Q), (img_R, img_S), (img_X, img_Y), retval)
}
