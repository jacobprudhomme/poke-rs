use core::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, projective_point::Point},
    polynomial_ring::poly::Polynomial,
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};

use crate::{
    SUCCESS_RETVAL,
    bn::BigNum,
    dlp::{solve_dlp_order_powers_of_two_three_five, solve_dlp_small_prime_power_order},
    endomorphism::{
        apply_endomorphism_from_quaternion, find_kernel_of_backtracking_isogeny_prime_power_degree,
        represent_integer,
    },
    masking::mask_basis_by_same_scalar,
    poke,
    rand::sample_random_torsion_basis,
};

pub fn generate_2d_isogeny_poke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_P: usize,
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
        q.widen() * &q_dual.widen() * &pub_params.three_times_five_torsion_order.widen();
    let (theta, ok) = represent_integer(&norm, &pub_params.field_characteristic);
    retval &= ok;

    /* Construct backtracking isogeny */

    // Find kernels of degree-3^b and degree-5^c backtracking isogenies of theta
    let K3 = find_kernel_of_backtracking_isogeny_prime_power_degree(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &pub_params.three_torsion_basis,
        3,
        pub_params.three_torsion_exp,
        &pub_params.reduced_three_torsion_order,
        &pub_params.three_torsion_order,
        &pub_params.three_adic_basis,
    );
    let K5 = find_kernel_of_backtracking_isogeny_prime_power_degree(
        &pub_params.field_characteristic,
        &pub_params.starting_curve,
        &theta,
        &pub_params.five_torsion_basis,
        5,
        pub_params.five_torsion_exp,
        &pub_params.reduced_five_torsion_order,
        &pub_params.five_torsion_order,
        &pub_params.five_adic_basis,
    );
    let K35 = pub_params.starting_curve.add(&K3, &K5);

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

    let codomain = pub_params
        .starting_curve
        .velu_composite_isogeny::<Polynomial<Fp2>>(
            &K35.to_pointx(),
            &[
                (3, pub_params.three_torsion_exp),
                (5, pub_params.five_torsion_exp),
            ],
            &mut theta_two_torsion_basis,
        );

    let theta_backtrack_two_torsion_basis_img = BasisX::from_slice(&theta_two_torsion_basis);
    let (theta_backtrack_P, theta_backtrack_Q) =
        codomain.lift_basis(&theta_backtrack_two_torsion_basis_img);

    // Get rid of the (3^b * 5^c)-factor introduced by the backtracking isogeny
    let (scaled_theta_backtrack_P, scaled_theta_backtrack_Q) = mask_basis_by_same_scalar(
        &codomain,
        &(theta_backtrack_P, theta_backtrack_Q),
        &pub_params
            .three_times_five_torsion_order
            .invert_mod(&pub_params.full_two_torsion_order),
    );
    // TODO: Test if removing the negation of these points still creates a valid isogeny
    let (scaled_P, scaled_Q) = mask_basis_by_same_scalar(&pub_params.starting_curve, &(-P, -Q), &q);

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
    // Used to make a choice of scalar factor later
    // FIXME: if we're computing this in addition to the proper pairing, would it not be better to just
    // compute the pairing from the power of e(U,V) directly, and fix the same power in both keygen and decryption?
    let eUV_power_q = eUV_AB.pow(&degree.to_le_bytes(), degree.nbits());
    let eUV_aux_is_eUV_power_q = eUV_aux.equals(&eUV_power_q);

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
    img_P.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

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
    img_Q.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

    (aux_curve, (img_P, img_Q), retval)
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
