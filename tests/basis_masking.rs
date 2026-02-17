#![allow(non_snake_case)]

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point};
use poke::{
    SUCCESS_RETVAL,
    bn::BigNum,
    params::{poke_i, poke_iii, poke_v},
    poke::PublicParams,
    rand::{
        sample_random_invertible_matrix_mod_special_prime_power, sample_random_torsion_basis,
        sample_random_unit_mod_special_prime_power,
    },
};
use rstest::rstest;

fn multiply_basis_by_scalars<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &(Point<Fp2>, Point<Fp2>),
    s1: &BigNum<NUM_WORDS>,
    s2: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let masked_P = curve.mul(&basis.0, &s1.to_le_bytes(), s1.nbits());
    let masked_Q = curve.mul(&basis.1, &s2.to_le_bytes(), s2.nbits());
    let masked_PQ = curve.sub(&masked_P, &masked_Q);

    BasisX::from_points(
        &masked_P.to_pointx(),
        &masked_Q.to_pointx(),
        &masked_PQ.to_pointx(),
    )
}

fn multiply_basis_by_scalar_matrix<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &(Point<Fp2>, Point<Fp2>),
    S: &[[BigNum<NUM_WORDS>; 2]; 2],
) -> BasisX<Fp2> {
    let masked_P = curve.add(
        &curve.mul(&basis.0, &S[0][0].to_le_bytes(), S[0][0].nbits()),
        &curve.mul(&basis.1, &S[0][1].to_le_bytes(), S[0][1].nbits()),
    );
    let masked_Q = curve.add(
        &curve.mul(&basis.0, &S[1][0].to_le_bytes(), S[1][0].nbits()),
        &curve.mul(&basis.1, &S[1][1].to_le_bytes(), S[1][1].nbits()),
    );
    let masked_PQ = curve.sub(&masked_P, &masked_Q);

    BasisX::from_points(
        &masked_P.to_pointx(),
        &masked_Q.to_pointx(),
        &masked_PQ.to_pointx(),
    )
}

fn multiply_xonly_basis_by_same_scalar_xmul<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    s: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let masked_xP = curve.xmul(&basis.P, &s.to_le_bytes(), s.nbits());
    let masked_xQ = curve.xmul(&basis.Q, &s.to_le_bytes(), s.nbits());
    let masked_xPQ = curve.xmul(&basis.PQ, &s.to_le_bytes(), s.nbits());

    BasisX::from_points(&masked_xP, &masked_xQ, &masked_xPQ)
}

// fn multiply_xonly_basis_by_same_scalar_projective_difference<Fp2: Fp2Trait>(
//     curve: &Curve<Fp2>,
//     basis: &BasisX<Fp2>,
//     s: &BigNum,
// ) -> BasisX<Fp2> {
//     let masked_xP = curve.xmul(&basis.P, s.as_le_bytes(), s.nbits());
//     let masked_xPQ = curve.xmul(&basis.PQ, s.as_le_bytes(), s.nbits());
//     let masked_xQ = curve.projective_difference(&masked_xP, &masked_xPQ);

//     BasisX::from_points(&masked_xP, &masked_xQ, &masked_xPQ)
// }

fn multiply_xonly_basis_by_scalars_flipped_Q_basis<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    s1: &BigNum<NUM_WORDS>,
    s2: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let xPpQ = Curve::<Fp2>::xdiff_add(&basis.P, &basis.Q, &basis.PQ);
    let basis_flipped_Q = BasisX::from_points(&basis.P, &basis.Q, &xPpQ);

    let masked_xP = curve.xmul(&basis.P, &s1.to_le_bytes(), s1.nbits());
    let masked_xQ = curve.xmul(&basis.Q, &s2.to_le_bytes(), s2.nbits());
    let masked_xPQ = curve.ladder_biscalar(
        &basis_flipped_Q,
        &s1.to_le_bytes(),
        &s2.to_le_bytes(),
        s1.nbits(),
        s2.nbits(),
    );

    BasisX::from_points(&masked_xP, &masked_xQ, &masked_xPQ)
}

fn multiply_xonly_basis_by_scalars_negate_second_scalar<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    s1: &BigNum<NUM_WORDS>,
    s2: &BigNum<NUM_WORDS>,
    modulus: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let s2_neg = modulus - s2;

    let masked_xP = curve.xmul(&basis.P, &s1.to_le_bytes(), s1.nbits());
    let masked_xQ = curve.xmul(&basis.Q, &s2.to_le_bytes(), s2.nbits());
    let masked_xPQ = curve.ladder_biscalar(
        &basis,
        &s1.to_le_bytes(),
        &s2_neg.to_le_bytes(),
        s1.nbits(),
        s2_neg.nbits(),
    );

    BasisX::from_points(&masked_xP, &masked_xQ, &masked_xPQ)
}

fn multiply_xonly_basis_by_scalars_using_invert_first_scalar<
    Fp2: Fp2Trait,
    const NUM_WORDS: usize,
>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    s1: &BigNum<NUM_WORDS>,
    s2: &BigNum<NUM_WORDS>,
    modulus: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    // SAFETY: we expect units as input, just as in the protocol
    let s1_inv_neg = modulus - s1.invert_mod_vartime(&modulus);
    let s1_inv_neg_s2 = s1_inv_neg * s2;

    let masked_xP = curve.xmul(&basis.P, &s1.to_le_bytes(), s1.nbits());
    let masked_xQ = curve.xmul(&basis.Q, &s2.to_le_bytes(), s2.nbits());

    let scaled_xQ = curve.xmul(&basis.Q, &s1.to_le_bytes(), s1.nbits());
    let scaled_xPQ = curve.xmul(&basis.PQ, &s1.to_le_bytes(), s1.nbits());
    let scaled_basis = BasisX::from_points(&masked_xP, &scaled_xQ, &scaled_xPQ);

    let masked_xPQ = curve.three_point_ladder(
        &scaled_basis,
        &s1_inv_neg_s2.to_le_bytes(),
        s1_inv_neg_s2.nbits(),
    );

    BasisX::from_points(&masked_xP, &masked_xQ, &masked_xPQ)
}

fn multiply_xonly_basis_by_scalar_matrix_scalar_difference<
    Fp2: Fp2Trait,
    const NUM_WORDS: usize,
>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    S: &[[BigNum<NUM_WORDS>; 2]; 2],
    modulus: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let s_diff1 = &S[0][0] + modulus - &S[1][0];
    let s_diff2 = &S[0][1] + modulus - &S[1][1];

    let masked_xP = curve.ladder_biscalar(
        basis,
        &S[0][0].to_le_bytes(),
        &S[0][1].to_le_bytes(),
        S[0][0].nbits(),
        S[0][1].nbits(),
    );
    let masked_xQ = curve.ladder_biscalar(
        basis,
        &S[1][0].to_le_bytes(),
        &S[1][1].to_le_bytes(),
        S[1][0].nbits(),
        S[1][1].nbits(),
    );
    let masked_xPQ = curve.ladder_biscalar(
        basis,
        &s_diff1.to_le_bytes(),
        &s_diff2.to_le_bytes(),
        s_diff1.nbits(),
        s_diff2.nbits(),
    );

    BasisX::from_points(&masked_xP, &masked_xQ, &masked_xPQ)
}

fn special_case_multiply_basis_by_scalar_matrix<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &(Point<Fp2>, Point<Fp2>),
    S: &[[BigNum<NUM_WORDS>; 2]; 2],
) -> (PointX<Fp2>, PointX<Fp2>) {
    let masked_P = curve.add(
        &curve.mul(&basis.0, &S[0][0].to_le_bytes(), S[0][0].nbits()),
        &curve.mul(&basis.1, &S[0][1].to_le_bytes(), S[0][1].nbits()),
    );
    let masked_Q = curve.add(
        &curve.mul(&basis.0, &S[1][0].to_le_bytes(), S[1][0].nbits()),
        &curve.mul(&basis.1, &S[1][1].to_le_bytes(), S[1][1].nbits()),
    );

    (masked_P.to_pointx(), masked_Q.to_pointx())
}

fn special_case_multiply_xonly_basis_by_scalar_matrix<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    S: &[[BigNum<NUM_WORDS>; 2]; 2],
) -> (PointX<Fp2>, PointX<Fp2>) {
    let masked_xP = curve.ladder_biscalar(
        basis,
        &S[0][0].to_le_bytes(),
        &S[0][1].to_le_bytes(),
        S[0][0].nbits(),
        S[0][1].nbits(),
    );
    let masked_xQ = curve.ladder_biscalar(
        basis,
        &S[1][0].to_le_bytes(),
        &S[1][1].to_le_bytes(),
        S[1][0].nbits(),
        S[1][1].nbits(),
    );

    (masked_xP, masked_xQ)
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn all_methods_for_single_scalar_are_equal<
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
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    #[case] params: PublicParams<
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
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
) {
    let s = sample_random_unit_mod_special_prime_power(
        params.five_torsion.base,
        &params.five_torsion.order,
    );
    let (P, Q, _) = sample_random_torsion_basis(
        &params.starting_curve,
        &[params.five_torsion.base],
        &[&params.five_torsion.reduced_order],
        &params.five_torsion.order,
        &params.five_torsion.cofactor,
    );
    let xP = P.to_pointx();
    let xQ = Q.to_pointx();
    let xPQ = params.starting_curve.sub(&P, &Q).to_pointx();

    let res_lift = {
        let basisx = BasisX::from_points(&xP, &xQ, &xPQ);
        let basis = params.starting_curve.lift_basis(&basisx);
        multiply_basis_by_scalars(&params.starting_curve, &basis, &s, &s)
    };
    let res_lift_normalized = {
        let basis = params
            .starting_curve
            .lift_basis_normalised(&xP.x(), &xQ.x(), &xPQ.x());
        multiply_basis_by_scalars(&params.starting_curve, &basis, &s, &s)
    };
    let res_xmul = {
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        multiply_xonly_basis_by_same_scalar_xmul(&params.starting_curve, &basis, &s)
    };
    // let res_projective_difference = {
    //     let basis = BasisX::from_points(&xP, &xQ, &xPQ);
    //     multiply_xonly_basis_by_same_scalar_projective_difference(
    //         &params.starting_curve,
    //         &basis,
    //         &s,
    //     )
    // };

    assert_eq!(res_lift.P.equals(&res_lift_normalized.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_lift_normalized.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_lift_normalized.PQ), SUCCESS_RETVAL);

    assert_eq!(res_lift.P.equals(&res_xmul.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_xmul.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_xmul.PQ), SUCCESS_RETVAL);

    // assert_eq!(
    //     res_lift.P.equals(&res_projective_difference.P),
    //     SUCCESS_RETVAL,
    // );
    // assert_eq!(
    //     res_lift.Q.equals(&res_projective_difference.Q),
    //     SUCCESS_RETVAL,
    // );
    // assert_eq!(
    //     res_lift.PQ.equals(&res_projective_difference.PQ),
    //     SUCCESS_RETVAL,
    // );
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn all_methods_for_different_scalars_are_equal<
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
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    #[case] params: PublicParams<
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
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
) {
    let s1 = sample_random_unit_mod_special_prime_power(
        params.five_torsion.base,
        &params.five_torsion.order,
    );
    let s2 = sample_random_unit_mod_special_prime_power(
        params.five_torsion.base,
        &params.five_torsion.order,
    );
    let (P, Q, _) = sample_random_torsion_basis(
        &params.starting_curve,
        &[params.five_torsion.base],
        &[&params.five_torsion.reduced_order],
        &params.five_torsion.order,
        &params.five_torsion.cofactor,
    );
    let xP = P.to_pointx();
    let xQ = Q.to_pointx();
    let xPQ = params.starting_curve.sub(&P, &Q).to_pointx();

    let res_lift = {
        let basisx = BasisX::from_points(&xP, &xQ, &xPQ);
        let basis = params.starting_curve.lift_basis(&basisx);
        multiply_basis_by_scalars(&params.starting_curve, &basis, &s1, &s2)
    };
    let res_lift_normalized = {
        let basis = params
            .starting_curve
            .lift_basis_normalised(&xP.x(), &xQ.x(), &xPQ.x());
        multiply_basis_by_scalars(&params.starting_curve, &basis, &s1, &s2)
    };
    let res_flipped_Q = {
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        multiply_xonly_basis_by_scalars_flipped_Q_basis(&params.starting_curve, &basis, &s1, &s2)
    };
    let res_neg_s2 = {
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        multiply_xonly_basis_by_scalars_negate_second_scalar(
            &params.starting_curve,
            &basis,
            &s1,
            &s2,
            &params.five_torsion.order,
        )
    };
    let res_inv_s1 = {
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        multiply_xonly_basis_by_scalars_using_invert_first_scalar(
            &params.starting_curve,
            &basis,
            &s1,
            &s2,
            &params.five_torsion.order,
        )
    };

    assert_eq!(res_lift.P.equals(&res_lift_normalized.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_lift_normalized.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_lift_normalized.PQ), SUCCESS_RETVAL);

    assert_eq!(res_lift.P.equals(&res_flipped_Q.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_flipped_Q.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_flipped_Q.PQ), SUCCESS_RETVAL);

    assert_eq!(res_lift.P.equals(&res_neg_s2.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_neg_s2.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_neg_s2.PQ), SUCCESS_RETVAL);

    assert_eq!(res_lift.P.equals(&res_inv_s1.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_inv_s1.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_inv_s1.PQ), SUCCESS_RETVAL);
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn all_methods_for_scalar_matrix_are_equal<
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
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    #[case] params: PublicParams<
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
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
) {
    let S = sample_random_invertible_matrix_mod_special_prime_power(
        params.five_torsion.base,
        &params.five_torsion.order,
    );
    let (P, Q, _) = sample_random_torsion_basis(
        &params.starting_curve,
        &[params.five_torsion.base],
        &[&params.five_torsion.reduced_order],
        &params.five_torsion.order,
        &params.five_torsion.cofactor,
    );
    let xP = P.to_pointx();
    let xQ = Q.to_pointx();
    let xPQ = params.starting_curve.sub(&P, &Q).to_pointx();

    let res_lift = {
        let basisx = BasisX::from_points(&xP, &xQ, &xPQ);
        let basis = params.starting_curve.lift_basis(&basisx);
        multiply_basis_by_scalar_matrix(&params.starting_curve, &basis, &S)
    };
    let res_lift_normalized = {
        let basis = params
            .starting_curve
            .lift_basis_normalised(&xP.x(), &xQ.x(), &xPQ.x());
        multiply_basis_by_scalar_matrix(&params.starting_curve, &basis, &S)
    };
    let res_scalar_difference = {
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        multiply_xonly_basis_by_scalar_matrix_scalar_difference(
            &params.starting_curve,
            &basis,
            &S,
            &params.five_torsion.order,
        )
    };

    assert_eq!(res_lift.P.equals(&res_lift_normalized.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_lift_normalized.Q), SUCCESS_RETVAL);
    assert_eq!(res_lift.PQ.equals(&res_lift_normalized.PQ), SUCCESS_RETVAL);

    assert_eq!(res_lift.P.equals(&res_scalar_difference.P), SUCCESS_RETVAL);
    assert_eq!(res_lift.Q.equals(&res_scalar_difference.Q), SUCCESS_RETVAL);
    assert_eq!(
        res_lift.PQ.equals(&res_scalar_difference.PQ),
        SUCCESS_RETVAL,
    );
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn all_methods_for_special_case_are_equal<
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
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    #[case] params: PublicParams<
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
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
) {
    let S = sample_random_invertible_matrix_mod_special_prime_power(
        params.five_torsion.base,
        &params.five_torsion.order,
    );
    let (P, Q, _) = sample_random_torsion_basis(
        &params.starting_curve,
        &[params.five_torsion.base],
        &[&params.five_torsion.reduced_order],
        &params.five_torsion.order,
        &params.five_torsion.cofactor,
    );
    let xP = P.to_pointx();
    let xQ = Q.to_pointx();
    let xPQ = params.starting_curve.sub(&P, &Q).to_pointx();

    let res_lift = {
        let basisx = BasisX::from_points(&xP, &xQ, &xPQ);
        let basis = params.starting_curve.lift_basis(&basisx);
        special_case_multiply_basis_by_scalar_matrix(&params.starting_curve, &basis, &S)
    };
    let res_lift_normalized = {
        let basis = params
            .starting_curve
            .lift_basis_normalised(&xP.x(), &xQ.x(), &xPQ.x());
        special_case_multiply_basis_by_scalar_matrix(&params.starting_curve, &basis, &S)
    };
    let res_biscalar_ladder = {
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        special_case_multiply_xonly_basis_by_scalar_matrix(&params.starting_curve, &basis, &S)
    };

    assert_eq!(res_lift.0.equals(&res_lift_normalized.0), SUCCESS_RETVAL);
    assert_eq!(res_lift.1.equals(&res_lift_normalized.1), SUCCESS_RETVAL);

    assert_eq!(res_lift.0.equals(&res_biscalar_ladder.0), SUCCESS_RETVAL);
    assert_eq!(res_lift.1.equals(&res_biscalar_ladder.1), SUCCESS_RETVAL);
}
