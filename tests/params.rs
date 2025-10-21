#![allow(non_snake_case)]

use fp2::traits::Fp;
use num_bigint::BigUint;
use poke::{
    PublicParams,
    fields::{PokeFieldI, PokeFieldIII, PokeFieldV},
    params::{
        poke_i::create_poke_i_params, poke_iii::create_poke_iii_params,
        poke_v::create_poke_v_params,
    },
};
use rstest::{fixture, rstest};

const SUCCESS_RETVAL: u32 = 0xffffffff;

#[fixture]
fn params_poke_i() -> PublicParams<PokeFieldI> {
    create_poke_i_params()
}

#[fixture]
fn params_poke_iii() -> PublicParams<PokeFieldIII> {
    create_poke_iii_params()
}

#[fixture]
fn params_poke_v() -> PublicParams<PokeFieldV> {
    create_poke_v_params()
}

#[rstest]
fn starting_curve_has_j_invariant_1728_poke_i(params_poke_i: PublicParams<PokeFieldI>) {
    assert_eq!(
        params_poke_i
            .starting_curve
            .j_invariant()
            .equals(&PokeFieldI::from_i32(1728)),
        SUCCESS_RETVAL,
        "j-invariant of E_0 is not 1728",
    );
}

#[rstest]
fn starting_curve_has_j_invariant_1728_poke_iii(params_poke_iii: PublicParams<PokeFieldIII>) {
    assert_eq!(
        params_poke_iii
            .starting_curve
            .j_invariant()
            .equals(&PokeFieldIII::from_i32(1728)),
        SUCCESS_RETVAL,
        "j-invariant of E_0 is not 1728",
    );
}

#[rstest]
fn starting_curve_has_j_invariant_1728_poke_v(params_poke_v: PublicParams<PokeFieldV>) {
    assert_eq!(
        params_poke_v
            .starting_curve
            .j_invariant()
            .equals(&PokeFieldV::from_i32(1728)),
        SUCCESS_RETVAL,
        "j-invariant of E_0 is not 1728",
    );
}

#[rstest]
fn torsion_basis_points_have_correct_order_poke_i(params_poke_i: PublicParams<PokeFieldI>) {
    // Represent 3^b and 5^c as byte arrays
    let three_torsion_order = BigUint::from(3u8).pow(
        params_poke_i
            .three_torsion_exp
            .try_into()
            .expect("Exponent of the 3-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let three_torsion_order_bitsize = three_torsion_order.bits().try_into().expect(
        "Size in bits of 3^b is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let three_torsion_order = three_torsion_order.to_bytes_le();

    let five_torsion_order = BigUint::from(5u8).pow(
        params_poke_i
            .five_torsion_exp
            .try_into()
            .expect("Exponent of the 5-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let five_torsion_order_bitsize = five_torsion_order.bits().try_into().expect(
        "Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let five_torsion_order = five_torsion_order.to_bytes_le();

    // Compute [2^a] * (P, Q, P - Q)
    let two_torsion_basis_times_order = params_poke_i.starting_curve.basis_double_iter(
        &params_poke_i.two_torsion_basis,
        params_poke_i.two_torsion_exp,
    );
    assert_eq!(two_torsion_basis_times_order.P.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.Q.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.PQ.is_zero(), SUCCESS_RETVAL);

    // Compute [3^b] * (R, S, R - S)
    let xR_times_order = params_poke_i.starting_curve.xmul(
        &params_poke_i.three_torsion_basis.P,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xS_times_order = params_poke_i.starting_curve.xmul(
        &params_poke_i.three_torsion_basis.Q,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xRS_times_order = params_poke_i.starting_curve.xmul(
        &params_poke_i.three_torsion_basis.PQ,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    assert_eq!(xR_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xS_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xRS_times_order.is_zero(), SUCCESS_RETVAL);

    // Compute [5^c] * (X, Y, X - Y)
    let xX_times_order = params_poke_i.starting_curve.xmul(
        &params_poke_i.five_torsion_basis.P,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xY_times_order = params_poke_i.starting_curve.xmul(
        &params_poke_i.five_torsion_basis.Q,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xXY_times_order = params_poke_i.starting_curve.xmul(
        &params_poke_i.five_torsion_basis.PQ,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    assert_eq!(xX_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xY_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xXY_times_order.is_zero(), SUCCESS_RETVAL);
}

#[rstest]
fn torsion_basis_points_have_correct_order_poke_iii(params_poke_iii: PublicParams<PokeFieldIII>) {
    // Represent 3^b and 5^c as byte arrays
    let three_torsion_order = BigUint::from(3u8).pow(
        params_poke_iii
            .three_torsion_exp
            .try_into()
            .expect("Exponent of the 3-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let three_torsion_order_bitsize = three_torsion_order.bits().try_into().expect(
        "Size in bits of 3^b is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let three_torsion_order = three_torsion_order.to_bytes_le();

    let five_torsion_order = BigUint::from(5u8).pow(
        params_poke_iii
            .five_torsion_exp
            .try_into()
            .expect("Exponent of the 5-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let five_torsion_order_bitsize = five_torsion_order.bits().try_into().expect(
        "Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let five_torsion_order = five_torsion_order.to_bytes_le();

    // Compute [2^a] * (P, Q, P - Q)
    let two_torsion_basis_times_order = params_poke_iii.starting_curve.basis_double_iter(
        &params_poke_iii.two_torsion_basis,
        params_poke_iii.two_torsion_exp,
    );
    assert_eq!(two_torsion_basis_times_order.P.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.Q.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.PQ.is_zero(), SUCCESS_RETVAL);

    // Compute [3^b] * (R, S, R - S)
    let xR_times_order = params_poke_iii.starting_curve.xmul(
        &params_poke_iii.three_torsion_basis.P,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xS_times_order = params_poke_iii.starting_curve.xmul(
        &params_poke_iii.three_torsion_basis.Q,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xRS_times_order = params_poke_iii.starting_curve.xmul(
        &params_poke_iii.three_torsion_basis.PQ,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    assert_eq!(xR_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xS_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xRS_times_order.is_zero(), SUCCESS_RETVAL);

    // Compute [5^c] * (X, Y, X - Y)
    let xX_times_order = params_poke_iii.starting_curve.xmul(
        &params_poke_iii.five_torsion_basis.P,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xY_times_order = params_poke_iii.starting_curve.xmul(
        &params_poke_iii.five_torsion_basis.Q,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xXY_times_order = params_poke_iii.starting_curve.xmul(
        &params_poke_iii.five_torsion_basis.PQ,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    assert_eq!(xX_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xY_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xXY_times_order.is_zero(), SUCCESS_RETVAL);
}

#[rstest]
fn torsion_basis_points_have_correct_order_poke_v(params_poke_v: PublicParams<PokeFieldV>) {
    // Represent 3^b and 5^c as byte arrays
    let three_torsion_order = BigUint::from(3u8).pow(
        params_poke_v
            .three_torsion_exp
            .try_into()
            .expect("Exponent of the 3-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let three_torsion_order_bitsize = three_torsion_order.bits().try_into().expect(
        "Size in bits of 3^b is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let three_torsion_order = three_torsion_order.to_bytes_le();

    let five_torsion_order = BigUint::from(5u8).pow(
        params_poke_v
            .five_torsion_exp
            .try_into()
            .expect("Exponent of the 5-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let five_torsion_order_bitsize = five_torsion_order.bits().try_into().expect(
        "Size in bits of 5^c is too big to fit in a usize (we do not ever expect this to happen)",
    );
    let five_torsion_order = five_torsion_order.to_bytes_le();

    // Compute [2^a] * (P, Q, P - Q)
    let two_torsion_basis_times_order = params_poke_v.starting_curve.basis_double_iter(
        &params_poke_v.two_torsion_basis,
        params_poke_v.two_torsion_exp,
    );
    assert_eq!(two_torsion_basis_times_order.P.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.Q.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.PQ.is_zero(), SUCCESS_RETVAL);

    // Compute [3^b] * (R, S, R - S)
    let xR_times_order = params_poke_v.starting_curve.xmul(
        &params_poke_v.three_torsion_basis.P,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xS_times_order = params_poke_v.starting_curve.xmul(
        &params_poke_v.three_torsion_basis.Q,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    let xRS_times_order = params_poke_v.starting_curve.xmul(
        &params_poke_v.three_torsion_basis.PQ,
        &three_torsion_order,
        three_torsion_order_bitsize,
    );
    assert_eq!(xR_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xS_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xRS_times_order.is_zero(), SUCCESS_RETVAL);

    // Compute [5^c] * (X, Y, X - Y)
    let xX_times_order = params_poke_v.starting_curve.xmul(
        &params_poke_v.five_torsion_basis.P,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xY_times_order = params_poke_v.starting_curve.xmul(
        &params_poke_v.five_torsion_basis.Q,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    let xXY_times_order = params_poke_v.starting_curve.xmul(
        &params_poke_v.five_torsion_basis.PQ,
        &five_torsion_order,
        five_torsion_order_bitsize,
    );
    assert_eq!(xX_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xY_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xXY_times_order.is_zero(), SUCCESS_RETVAL);
}

#[rstest]
fn torsion_basis_points_are_on_curve_poke_i(params_poke_i: PublicParams<PokeFieldI>) {
    // Check (P, Q, P - Q), a basis of the 2^a-torsion
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.two_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "P is not on the curve E_0",
    );
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.two_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Q is not on the curve E_0",
    );
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.two_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "P - Q is not on the curve E_0",
    );

    // Check (R, S, R - S), a basis of the 3^b-torsion
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.three_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.three_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.three_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "R - S is not on the curve E_0",
    );

    // Check (X, Y, X - Y), a basis of the 5^c-torsion
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.five_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "X is not on the curve E_0",
    );
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.five_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Y is not on the curve E_0",
    );
    assert_eq!(
        params_poke_i
            .starting_curve
            .is_on_curve(&params_poke_i.five_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "X - Y is not on the curve E_0",
    );
}

#[rstest]
fn torsion_basis_points_are_on_curve_poke_iii(params_poke_iii: PublicParams<PokeFieldIII>) {
    // Check (P, Q, P - Q), a basis of the 2^a-torsion
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.two_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "P is not on the curve E_0",
    );
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.two_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Q is not on the curve E_0",
    );
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.two_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "P - Q is not on the curve E_0",
    );

    // Check (R, S, R - S), a basis of the 3^b-torsion
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.three_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.three_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.three_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "R - S is not on the curve E_0",
    );

    // Check (X, Y, X - Y), a basis of the 5^c-torsion
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.five_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "X is not on the curve E_0",
    );
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.five_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Y is not on the curve E_0",
    );
    assert_eq!(
        params_poke_iii
            .starting_curve
            .is_on_curve(&params_poke_iii.five_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "X - Y is not on the curve E_0",
    );
}

#[rstest]
fn torsion_basis_points_are_on_curve_poke_v(params_poke_v: PublicParams<PokeFieldV>) {
    // Check (P, Q, P - Q), a basis of the 2^a-torsion
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.two_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "P is not on the curve E_0",
    );
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.two_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Q is not on the curve E_0",
    );
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.two_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "P - Q is not on the curve E_0",
    );

    // Check (R, S, R - S), a basis of the 3^b-torsion
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.three_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.three_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "R is not on the curve E_0",
    );
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.three_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "R - S is not on the curve E_0",
    );

    // Check (X, Y, X - Y), a basis of the 5^c-torsion
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.five_torsion_basis.P.x()),
        SUCCESS_RETVAL,
        "X is not on the curve E_0",
    );
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.five_torsion_basis.Q.x()),
        SUCCESS_RETVAL,
        "Y is not on the curve E_0",
    );
    assert_eq!(
        params_poke_v
            .starting_curve
            .is_on_curve(&params_poke_v.five_torsion_basis.PQ.x()),
        SUCCESS_RETVAL,
        "X - Y is not on the curve E_0",
    );
}
