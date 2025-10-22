#![allow(non_snake_case)]

use fp2::traits::Fp2 as FpTrait;
use num_bigint::BigUint;
use poke::{
    PublicParams, SUCCESS_RETVAL,
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
    // Represent 3^b and 5^c as byte arrays
    let three_torsion_order = BigUint::from(3u8).pow(
        params
            .three_torsion_exp
            .try_into()
            .expect("Exponent of the 3-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let three_torsion_order_bitsize = three_torsion_order.bits().try_into().expect(
        "Size in bits of 3^b is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let three_torsion_order = three_torsion_order.to_bytes_le();

    let five_torsion_order = BigUint::from(5u8).pow(
        params
            .five_torsion_exp
            .try_into()
            .expect("Exponent of the 5-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let five_torsion_order_bitsize = five_torsion_order.bits().try_into().expect(
        "Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let five_torsion_order = five_torsion_order.to_bytes_le();

    // Compute [2^a] * (P, Q, P - Q)
    let two_torsion_basis_times_order = params
        .starting_curve
        .basis_double_iter(&params.two_torsion_basis, params.two_torsion_exp);
    assert_eq!(two_torsion_basis_times_order.P.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.Q.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.PQ.is_zero(), SUCCESS_RETVAL);

    // Compute [3^b] * (R, S, R - S)
    let xR_times_order = params.starting_curve.xmul(
        &params.three_torsion_basis.P,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xS_times_order = params.starting_curve.xmul(
        &params.three_torsion_basis.Q,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xRS_times_order = params.starting_curve.xmul(
        &params.three_torsion_basis.PQ,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    assert_eq!(xR_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xS_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xRS_times_order.is_zero(), SUCCESS_RETVAL);

    // Compute [5^c] * (X, Y, X - Y)
    let xX_times_order = params.starting_curve.xmul(
        &params.five_torsion_basis.P,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xY_times_order = params.starting_curve.xmul(
        &params.five_torsion_basis.Q,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xXY_times_order = params.starting_curve.xmul(
        &params.five_torsion_basis.PQ,
        &five_torsion_order,
        five_torsion_order_bitsize,
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
        "R is not on the curve E_0",
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
