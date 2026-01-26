#![allow(non_snake_case)]

use core::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};
use poke::{
    FAILURE_RETVAL, SUCCESS_RETVAL,
    bn::BigNum,
    dimtwo::{
        eval_2d_two_isogeny_chain_inke, eval_2d_two_isogeny_chain_inke_separate_bases,
        eval_2d_two_isogeny_chain_poke, eval_2d_two_isogeny_chain_poke_separate_bases,
        generate_2d_isogeny_inke, generate_2d_isogeny_poke,
    },
    inke::PublicParams as PublicParamsInke,
    masking::mask_basis_by_same_scalar,
    params,
    poke::PublicParams as PublicParamsPoke,
    rand::sample_random_secret_degree,
};
use rstest::rstest;

fn apply_2d_isogeny_to_full_torsion_basis_inke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_233: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
>(
    pub_params: &PublicParamsInke<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_233,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
    >,
    domain: &EllipticProduct<Fp2>,
    kernel: &(ProductPoint<Fp2>, ProductPoint<Fp2>),
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
) -> (
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
) {
    let codomain_curve = domain.curves().1;
    let (P1P2, Q1Q2) = kernel;

    let (P, Q) = pub_params
        .starting_curve
        .lift_basis(&pub_params.two_torsion_basis);
    let (R, S) = pub_params
        .starting_curve
        .lift_basis(&pub_params.three_torsion_basis);

    let mut PR = Point::INFINITY;
    pub_params.starting_curve.addto(&mut PR, &P);
    pub_params.starting_curve.addto(&mut PR, &R);
    let mut QS = Point::INFINITY;
    pub_params.starting_curve.addto(&mut QS, &Q);
    pub_params.starting_curve.addto(&mut QS, &S);
    let PRQS = pub_params.starting_curve.sub(&PR, &QS);

    let full_torsion_basis =
        BasisX::from_points(&PR.to_pointx(), &QS.to_pointx(), &PRQS.to_pointx());

    let (full_torsion_basis_EA, intermediate_curve, full_torsion_basis_EA1, _) =
        eval_2d_two_isogeny_chain_inke(
            domain,
            (P1P2, Q1Q2),
            pub_params.effective_two_torsion_exp,
            degree,
            degree_dual,
            (&pub_params.two_torsion, &pub_params.three_torsion),
            &full_torsion_basis,
            &pub_params.full_torsion_order,
            &pub_params.cofactor,
            PhantomData::<([(); NUM_WORDS_223], [(); NUM_WORDS_233])>,
        );

    let two_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &full_torsion_basis_EA,
        &pub_params.two_torsion.coproduct,
    );
    let two_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &two_torsion_basis_EA,
        &pub_params.two_torsion.inv_coproduct,
    );

    let three_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &full_torsion_basis_EA,
        &pub_params.three_torsion.coproduct,
    );
    let three_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &three_torsion_basis_EA,
        &pub_params.three_torsion.inv_coproduct,
    );

    let three_torsion_basis_EA1 = mask_basis_by_same_scalar(
        &intermediate_curve,
        &full_torsion_basis_EA1,
        &pub_params.three_torsion.coproduct,
    );
    let three_torsion_basis_EA1 = mask_basis_by_same_scalar(
        &intermediate_curve,
        &three_torsion_basis_EA1,
        &pub_params.three_torsion.inv_coproduct,
    );

    (
        two_torsion_basis_EA,
        three_torsion_basis_EA,
        three_torsion_basis_EA1,
    )
}

fn apply_2d_isogeny_to_individual_torsion_bases_inke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_233: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
>(
    pub_params: &PublicParamsInke<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_233,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
    >,
    domain: &EllipticProduct<Fp2>,
    kernel: &(ProductPoint<Fp2>, ProductPoint<Fp2>),
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
) -> (
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
) {
    let (P1P2, Q1Q2) = kernel;

    let (two_torsion_basis_EA, three_torsion_basis_EA, _, three_torsion_basis_EA1, _) =
        eval_2d_two_isogeny_chain_inke_separate_bases(
            domain,
            (P1P2, Q1Q2),
            pub_params.effective_two_torsion_exp,
            degree,
            degree_dual,
            (&pub_params.two_torsion, &pub_params.three_torsion),
            (
                &pub_params.two_torsion_basis,
                &pub_params.three_torsion_basis,
            ),
        );

    (
        two_torsion_basis_EA,
        three_torsion_basis_EA,
        three_torsion_basis_EA1,
    )
}

fn apply_2d_isogeny_to_full_torsion_basis_poke<
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
    pub_params: &PublicParamsPoke<
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
    domain: &EllipticProduct<Fp2>,
    kernel: &(ProductPoint<Fp2>, ProductPoint<Fp2>),
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
) -> (
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
) {
    let codomain_curve = domain.curves().1;
    let (P1P2, Q1Q2) = kernel;

    let (P, Q) = pub_params
        .starting_curve
        .lift_basis(&pub_params.two_torsion_basis);
    let (R, S) = pub_params
        .starting_curve
        .lift_basis(&pub_params.three_torsion_basis);
    let (X, Y) = pub_params
        .starting_curve
        .lift_basis(&pub_params.five_torsion_basis);

    let mut PRX = Point::INFINITY;
    pub_params.starting_curve.addto(&mut PRX, &P);
    pub_params.starting_curve.addto(&mut PRX, &R);
    pub_params.starting_curve.addto(&mut PRX, &X);
    let mut QSY = Point::INFINITY;
    pub_params.starting_curve.addto(&mut QSY, &Q);
    pub_params.starting_curve.addto(&mut QSY, &S);
    pub_params.starting_curve.addto(&mut QSY, &Y);
    let PRXQSY = pub_params.starting_curve.sub(&PRX, &QSY);

    let full_torsion_basis =
        BasisX::from_points(&PRX.to_pointx(), &QSY.to_pointx(), &PRXQSY.to_pointx());

    let (full_torsion_basis_EA, _) = eval_2d_two_isogeny_chain_poke(
        domain,
        (P1P2, Q1Q2),
        pub_params.effective_two_torsion_exp,
        degree,
        degree_dual,
        (
            &pub_params.two_torsion,
            &pub_params.three_torsion,
            &pub_params.five_torsion,
        ),
        &full_torsion_basis,
        &pub_params.full_torsion_order,
        &pub_params.cofactor,
        PhantomData::<(
            [(); NUM_WORDS_2235],
            [(); NUM_WORDS_2335],
            [(); NUM_WORDS_2355],
        )>,
    );

    let two_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &full_torsion_basis_EA,
        &pub_params.two_torsion.coproduct,
    );
    let two_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &two_torsion_basis_EA,
        &pub_params.two_torsion.inv_coproduct,
    );

    let three_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &full_torsion_basis_EA,
        &pub_params.three_torsion.coproduct,
    );
    let three_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &three_torsion_basis_EA,
        &pub_params.three_torsion.inv_coproduct,
    );

    let five_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &full_torsion_basis_EA,
        &pub_params.five_torsion.coproduct,
    );
    let five_torsion_basis_EA = mask_basis_by_same_scalar(
        &codomain_curve,
        &five_torsion_basis_EA,
        &pub_params.five_torsion.inv_coproduct,
    );

    (
        two_torsion_basis_EA,
        three_torsion_basis_EA,
        five_torsion_basis_EA,
    )
}

fn apply_2d_isogeny_to_individual_torsion_bases_poke<
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
    pub_params: &PublicParamsPoke<
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
    domain: &EllipticProduct<Fp2>,
    kernel: &(ProductPoint<Fp2>, ProductPoint<Fp2>),
    degree: &BigNum<NUM_WORDS_2>,
    degree_dual: &BigNum<NUM_WORDS_2>,
) -> (
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
    (Point<Fp2>, Point<Fp2>),
) {
    let (P1P2, Q1Q2) = kernel;

    let (two_torsion_basis_EA, three_torsion_basis_EA, five_torsion_basis_EA, _) =
        eval_2d_two_isogeny_chain_poke_separate_bases(
            domain,
            (P1P2, Q1Q2),
            pub_params.effective_two_torsion_exp,
            degree,
            degree_dual,
            (
                &pub_params.two_torsion,
                &pub_params.three_torsion,
                &pub_params.five_torsion,
            ),
            (
                &pub_params.two_torsion_basis,
                &pub_params.three_torsion_basis,
                &pub_params.five_torsion_basis,
            ),
        );

    (
        two_torsion_basis_EA,
        three_torsion_basis_EA,
        five_torsion_basis_EA,
    )
}

#[rstest]
#[case::inke_level_i(params::inke_i::get_params())]
#[case::inke_level_iii(params::inke_iii::get_params())]
#[case::inke_level_v(params::inke_v::get_params())]
fn both_methods_give_same_result_inke<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_233: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
>(
    #[case] params: PublicParamsInke<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_233,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
    >,
) {
    let (q, q_dual) = sample_random_secret_degree(&params.effective_two_torsion_order, &[3]);
    let (mut domain, mut kernel, mut ok) = generate_2d_isogeny_inke(&params, &q, &q_dual);
    while ok == FAILURE_RETVAL {
        (domain, kernel, ok) = generate_2d_isogeny_inke(&params, &q, &q_dual);
    }
    let ((P1, Q1), (R1, S1), (X1, Y1)) =
        apply_2d_isogeny_to_full_torsion_basis_inke(&params, &domain, &kernel, &q, &q_dual);
    let ((P2, Q2), (R2, S2), (X2, Y2)) =
        apply_2d_isogeny_to_individual_torsion_bases_inke(&params, &domain, &kernel, &q, &q_dual);

    assert_eq!(P1.to_pointx().equals(&P2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(Q1.to_pointx().equals(&Q2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(R1.to_pointx().equals(&R2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(S1.to_pointx().equals(&S2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(X1.to_pointx().equals(&X2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(Y1.to_pointx().equals(&Y2.to_pointx()), SUCCESS_RETVAL);
}

#[rstest]
#[case::poke_level_i(params::poke_i::get_params())]
#[case::poke_level_iii(params::poke_iii::get_params())]
#[case::poke_level_v(params::poke_v::get_params())]
fn both_methods_give_same_result_poke<
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
    #[case] params: PublicParamsPoke<
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
    let (q, q_dual) = sample_random_secret_degree(&params.effective_two_torsion_order, &[3, 5]);
    let (mut domain, mut kernel, mut ok) = generate_2d_isogeny_poke(&params, &q, &q_dual);
    while ok == FAILURE_RETVAL {
        (domain, kernel, ok) = generate_2d_isogeny_poke(&params, &q, &q_dual);
    }
    let ((P1, Q1), (R1, S1), (X1, Y1)) =
        apply_2d_isogeny_to_full_torsion_basis_poke(&params, &domain, &kernel, &q, &q_dual);
    let ((P2, Q2), (R2, S2), (X2, Y2)) =
        apply_2d_isogeny_to_individual_torsion_bases_poke(&params, &domain, &kernel, &q, &q_dual);

    assert_eq!(P1.to_pointx().equals(&P2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(Q1.to_pointx().equals(&Q2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(R1.to_pointx().equals(&R2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(S1.to_pointx().equals(&S2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(X1.to_pointx().equals(&X2.to_pointx()), SUCCESS_RETVAL);
    assert_eq!(Y1.to_pointx().equals(&Y2.to_pointx()), SUCCESS_RETVAL);
}
