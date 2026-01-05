use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};

use crate::{
    SUCCESS_RETVAL, bn::BigNum, dlp::solve_dlp_small_prime_power_order,
    rand::sample_random_torsion_basis_order_product_of_powers_of_small_primes,
};

pub fn eval_2d_two_isogeny_chain_on_prime_power_torsion_basis<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_ORD: usize,
    const NUM_WORDS_COF: usize,
>(
    domain: &EllipticProduct<Fp2>,
    kernel: (&ProductPoint<Fp2>, &ProductPoint<Fp2>),
    chain_length: usize,
    degree: &BigNum<NUM_WORDS_2>, // 2^chain_length
    embedded_isogeny_degree: &BigNum<NUM_WORDS_2>,
    torsion_basis: &BasisX<Fp2>,
    torsion_basis_order_base: u8,
    torsion_basis_order_exp: usize,
    torsion_basis_order: &BigNum<NUM_WORDS_ORD>,
    torsion_basis_cofactor: &BigNum<NUM_WORDS_COF>,
    p_adic_basis_for_torsion_basis_order_base: &[BigNum<NUM_WORDS_ORD>],
) -> ((Point<Fp2>, Point<Fp2>), u32) {
    let mut retval = SUCCESS_RETVAL;

    let (embedded_isogeny_domain, embedded_isogeny_codomain) = domain.curves();
    let (P1P2, Q1Q2) = kernel;

    // Factor that shows up in the application of the 2D-isogeny, from the dual that appears in its representation
    let dual_factor = degree - embedded_isogeny_degree;

    // Lift basis to full points, as the 2D-isogeny function requires this as input
    let (P, Q) = embedded_isogeny_domain.lift_basis(torsion_basis);
    let PQ = embedded_isogeny_domain.sub(&P, &Q);

    // Generate random basis of the 5^c-torsion on E_AB
    let (U, V, eUV_AB) = sample_random_torsion_basis_order_product_of_powers_of_small_primes(
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
    let eUV_power_q = eUV_AB.pow(
        &embedded_isogeny_degree.to_le_bytes(),
        embedded_isogeny_degree.nbits(),
    );
    let eUV_aux_is_eUV_power_q = eUV_aux.equals(&eUV_power_q);

    let eXV = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &PV_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let eXmU = aux_curve.weil_pairing(
        &P_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &PmU_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let eYV = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &V_aux_curve.to_pointx().x(),
        &QV_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );
    let eYmU = aux_curve.weil_pairing(
        &Q_aux_curve.to_pointx().x(),
        &U_aux_curve.to_pointx().x(),
        &QmU_aux_curve.to_pointx().x(),
        &torsion_basis_order.to_le_bytes(),
        torsion_basis_order.nbits(),
    );

    // Solve discrete logarithm between pairings to obtain expression of P' in terms of <U',V'>
    let (x, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eXV,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;
    let (y, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eXmU,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;
    let (w, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eYV,
        torsion_basis_order_base,
        torsion_basis_order_exp,
        p_adic_basis_for_torsion_basis_order_base,
    );
    retval &= ok;
    let (z, ok) = solve_dlp_small_prime_power_order(
        &eUV_aux,
        &eYmU,
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
        &embedded_isogeny_degree.to_le_bytes(),
        embedded_isogeny_degree.nbits(),
    );
    // [2^(a-2) - q] * P'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &dual_factor.to_le_bytes(),
        dual_factor.nbits(),
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
        &embedded_isogeny_degree.to_le_bytes(),
        embedded_isogeny_degree.nbits(),
    );
    // [2^(a-2) - q] * Q'
    embedded_isogeny_codomain.mul_into(
        &mut V_aux_curve,
        &PQ_aux_curve,
        &dual_factor.to_le_bytes(),
        dual_factor.nbits(),
    );
    // Choose the appropriate one of the above points depending on whether
    // e(U',V') = e(U,V)^q or e(U',V') = e(U,V)^(2^(a-2) - q)
    let mut img_Q = U_aux_curve;
    img_Q.set_cond(&V_aux_curve, !eUV_aux_is_eUV_power_q);

    ((img_P, img_Q), retval)
}
