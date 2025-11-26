#![allow(non_snake_case)]

use fp2::traits::Fp2 as FpTrait;
use num_bigint::BigUint;
use poke::{
    FAILURE_RETVAL, PublicParams, SUCCESS_RETVAL,
    params::{poke_i, poke_iii, poke_v},
};
use rstest::rstest;

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn starting_curve_has_j_invariant_1728<Fp2: FpTrait>(#[case] params: PublicParams<Fp2>) {
    assert_eq!(
        params
            .starting_curve
            .j_invariant()
            .equals(&Fp2::from_i32(1728)),
        SUCCESS_RETVAL,
        "j-invariant of E_0 is not 1728",
    );
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn torsion_basis_points_have_correct_order<Fp2: FpTrait>(#[case] params: PublicParams<Fp2>) {
    let ONE = BigUint::from(1u8);
    let THREE = BigUint::from(3u8);
    let FIVE = BigUint::from(5u8);

    // Ensure [2^a] * (P, Q, P - Q) = O, and [2^a - 1] * (P, Q, P - Q) != O
    let two_torsion_basis_times_one_less_than_order = params
        .starting_curve
        .basis_double_iter(&params.two_torsion_basis, params.full_two_torsion_exp - 1);
    assert_eq!(
        two_torsion_basis_times_one_less_than_order.P.is_zero(),
        FAILURE_RETVAL
    );
    assert_eq!(
        two_torsion_basis_times_one_less_than_order.Q.is_zero(),
        FAILURE_RETVAL
    );
    assert_eq!(
        two_torsion_basis_times_one_less_than_order.PQ.is_zero(),
        FAILURE_RETVAL
    );

    let two_torsion_basis_times_order = params
        .starting_curve
        .basis_double_iter(&params.two_torsion_basis, params.full_two_torsion_exp);
    assert_eq!(two_torsion_basis_times_order.P.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.Q.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.PQ.is_zero(), SUCCESS_RETVAL);

    // Ensure [3^b] * (R, S, R - S) = O, and [3^b - 1] * (R, S, R - S) != O
    let xR_times_one_less_than_order = params.starting_curve.xmul(
        &params.three_torsion_basis.P,
        &(&params.three_torsion_order / &THREE).to_bytes_le(),
        (&params.three_torsion_order / &THREE).bits() as usize,
    );
    let xS_times_one_less_than_order = params.starting_curve.xmul(
        &params.three_torsion_basis.Q,
        &(&params.three_torsion_order / &THREE).to_bytes_le(),
        (&params.three_torsion_order / &THREE).bits() as usize,
    );
    let xRS_times_one_less_than_order = params.starting_curve.xmul(
        &params.three_torsion_basis.PQ,
        &(&params.three_torsion_order / &THREE).to_bytes_le(),
        (&params.three_torsion_order / &THREE).bits() as usize,
    );
    assert_eq!(xR_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xS_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xRS_times_one_less_than_order.is_zero(), FAILURE_RETVAL);

    let xR_times_order = params.starting_curve.xmul(
        &params.three_torsion_basis.P,
        &params.three_torsion_order.to_bytes_le(),
        params.three_torsion_order.bits() as usize,
    );
    let xS_times_order = params.starting_curve.xmul(
        &params.three_torsion_basis.Q,
        &params.three_torsion_order.to_bytes_le(),
        params.three_torsion_order.bits() as usize,
    );
    let xRS_times_order = params.starting_curve.xmul(
        &params.three_torsion_basis.PQ,
        &params.three_torsion_order.to_bytes_le(),
        params.three_torsion_order.bits() as usize,
    );
    assert_eq!(xR_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xS_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xRS_times_order.is_zero(), SUCCESS_RETVAL);

    // Ensure [5^c] * (X, Y, X - Y) = O and [5^c - 1] * (X, Y, X - Y) != O
    let xX_times_one_less_than_order = params.starting_curve.xmul(
        &params.five_torsion_basis.P,
        &(&params.five_torsion_order / &FIVE).to_bytes_le(),
        (&params.five_torsion_order / &FIVE).bits() as usize,
    );
    let xY_times_one_less_than_order = params.starting_curve.xmul(
        &params.five_torsion_basis.Q,
        &(&params.five_torsion_order / &FIVE).to_bytes_le(),
        (&params.five_torsion_order / &FIVE).bits() as usize,
    );
    let xXY_times_one_less_than_order = params.starting_curve.xmul(
        &params.five_torsion_basis.PQ,
        &(&params.five_torsion_order / &FIVE).to_bytes_le(),
        (&params.five_torsion_order / &FIVE).bits() as usize,
    );
    assert_eq!(xX_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xY_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xXY_times_one_less_than_order.is_zero(), FAILURE_RETVAL);

    let xX_times_order = params.starting_curve.xmul(
        &params.five_torsion_basis.P,
        &params.five_torsion_order.to_bytes_le(),
        params.five_torsion_order.bits() as usize,
    );
    let xY_times_order = params.starting_curve.xmul(
        &params.five_torsion_basis.Q,
        &params.five_torsion_order.to_bytes_le(),
        params.five_torsion_order.bits() as usize,
    );
    let xXY_times_order = params.starting_curve.xmul(
        &params.five_torsion_basis.PQ,
        &params.five_torsion_order.to_bytes_le(),
        params.five_torsion_order.bits() as usize,
    );
    assert_eq!(xX_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xY_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xXY_times_order.is_zero(), SUCCESS_RETVAL);
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn torsion_basis_points_are_on_curve<Fp2: FpTrait>(#[case] params: PublicParams<Fp2>) {
    // Check (P, Q, P - Q), a basis of the 2^a-torsion
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.two_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "P is not on the curve E_0",
    );
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.two_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Q is not on the curve E_0",
    );
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.two_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "P - Q is not on the curve E_0",
    );

    // Check (R, S, R - S), a basis of the 3^b-torsion
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.three_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.three_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "S is not on the curve E_0",
    );
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.three_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "R - S is not on the curve E_0",
    );

    // Check (X, Y, X - Y), a basis of the 5^c-torsion
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.five_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "X is not on the curve E_0",
    );
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.five_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Y is not on the curve E_0",
    );
    assert_eq!(
        params
            .starting_curve
            .is_on_curve(&params.five_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "X - Y is not on the curve E_0",
    );
}

#[rstest]
#[case::poke_level_i(poke_i::get_params())]
#[case::poke_level_iii(poke_iii::get_params())]
#[case::poke_level_v(poke_v::get_params())]
fn torsion_basis_points_are_linearly_independent<Fp2: FpTrait>(#[case] params: PublicParams<Fp2>) {
    let ONE = BigUint::from(1u8);
    let TWO = BigUint::from(2u8);
    let THREE = BigUint::from(3u8);
    let FIVE = BigUint::from(5u8);

    // Check (P, Q, P - Q), a basis of the 2^a-torsion
    let pair = params.starting_curve.weil_pairing_2exp(
        &params.two_torsion_basis.P.x(),
        &params.two_torsion_basis.Q.x(),
        &params.two_torsion_basis.PQ.x(),
        params.full_two_torsion_exp,
    );
    assert_eq!(
        pair.pow(
            &(&params.full_two_torsion_order / &TWO).to_bytes_le(),
            (&params.full_two_torsion_order / &TWO).bits() as usize
        )
        .equals(&Fp2::ONE),
        FAILURE_RETVAL,
    );
    assert_eq!(
        pair.pow(
            &params.full_two_torsion_order.to_bytes_le(),
            params.full_two_torsion_order.bits() as usize
        )
        .equals(&Fp2::ONE),
        SUCCESS_RETVAL,
    );

    // Check (R, S, R - S), a basis of the 3^b-torsion
    let pair = params.starting_curve.weil_pairing(
        &params.three_torsion_basis.P.x(),
        &params.three_torsion_basis.Q.x(),
        &params.three_torsion_basis.PQ.x(),
        &params.three_torsion_order.to_bytes_le(),
        params.three_torsion_order.bits() as usize,
    );
    assert_eq!(
        pair.pow(
            &(&params.three_torsion_order / &THREE).to_bytes_le(),
            (&params.three_torsion_order / &THREE).bits() as usize
        )
        .equals(&Fp2::ONE),
        FAILURE_RETVAL,
    );
    assert_eq!(
        pair.pow(
            &params.three_torsion_order.to_bytes_le(),
            params.three_torsion_order.bits() as usize
        )
        .equals(&Fp2::ONE),
        SUCCESS_RETVAL,
    );

    // Check (X, Y, X - Y), a basis of the 5^c-torsion
    let pair = params.starting_curve.weil_pairing(
        &params.five_torsion_basis.P.x(),
        &params.five_torsion_basis.Q.x(),
        &params.five_torsion_basis.PQ.x(),
        &params.five_torsion_order.to_bytes_le(),
        params.five_torsion_order.bits() as usize,
    );
    assert_eq!(
        pair.pow(
            &(&params.five_torsion_order / &FIVE).to_bytes_le(),
            (&params.five_torsion_order / &FIVE).bits() as usize
        )
        .equals(&Fp2::ONE),
        FAILURE_RETVAL,
    );
    assert_eq!(
        pair.pow(
            &params.five_torsion_order.to_bytes_le(),
            params.five_torsion_order.bits() as usize
        )
        .equals(&Fp2::ONE),
        SUCCESS_RETVAL,
    );
}
