#![allow(non_snake_case)]

use isogeny::elliptic::basis::BasisX;
use poke::{
    SUCCESS_RETVAL,
    bn::BigNum,
    fields::{PokeFieldI, PokeFieldIBase},
    params::poke_i,
    poke::PublicParams,
    rand::sample_random_unit_mod,
};
use rstest::{fixture, rstest};

#[fixture]
fn params() -> PublicParams<PokeFieldI> {
    poke_i::get_params()
}

#[fixture]
fn scalars(params: PublicParams<PokeFieldI>) -> (BigNum, BigNum) {
    sample_random_unit_mod(&params.full_two_torsion_order)
}

fn method1(
    params: &PublicParams<PokeFieldI>,
    (s, s_inv): (&BigNum, &BigNum),
) -> BasisX<PokeFieldI> {
    let (P, Q) = params.starting_curve.lift_basis(&params.two_torsion_basis);

    let P = params.starting_curve.mul(&P, s.as_le_bytes(), s.nbits());
    let Q = params
        .starting_curve
        .mul(&Q, s_inv.as_le_bytes(), s_inv.nbits());

    let PQ = params.starting_curve.sub(&P, &Q);

    BasisX::from_points(&P.to_pointx(), &Q.to_pointx(), &PQ.to_pointx())
}

fn method2(
    params: &PublicParams<PokeFieldI>,
    (s, s_inv): (&BigNum, &BigNum),
) -> BasisX<PokeFieldI> {
    let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();

    let P_x = params.starting_curve.xmul(&P_x, s.as_le_bytes(), s.nbits());
    let Q_x = params
        .starting_curve
        .xmul(&Q_x, s_inv.as_le_bytes(), s_inv.nbits());

    let (P, _) = params.starting_curve.lift_pointx(&P_x);
    let (Q, _) = params.starting_curve.lift_pointx(&Q_x);

    let PQ = params.starting_curve.sub(&P, &Q);

    BasisX::from_points(&P_x, &Q_x, &PQ.to_pointx())
}

fn method3(
    params: &PublicParams<PokeFieldI>,
    (s, s_inv): (&BigNum, &BigNum),
) -> BasisX<PokeFieldI> {
    let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();

    let P_x = params.starting_curve.xmul(&P_x, s.as_le_bytes(), s.nbits());
    let Q_x = params
        .starting_curve
        .xmul(&Q_x, s_inv.as_le_bytes(), s_inv.nbits());

    let (s_inv, _) = PokeFieldIBase::decode(s_inv.as_le_bytes());
    let minus_s_inv = (PokeFieldIBase::MINUS_ONE * s_inv).encode();

    let PQ_x = params.starting_curve.ladder_biscalar(
        &params.two_torsion_basis,
        s.as_le_bytes(),
        &minus_s_inv,
        s.nbits(),
        PokeFieldIBase::ENCODED_LENGTH,
    );

    BasisX::from_points(&P_x, &Q_x, &PQ_x)
}

fn projective_difference_method(
    params: &PublicParams<PokeFieldI>,
    (s, s_inv): (&BigNum, &BigNum),
) -> BasisX<PokeFieldI> {
    let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();

    let P_x = params.starting_curve.xmul(&P_x, s.as_le_bytes(), s.nbits());
    let Q_x = params
        .starting_curve
        .xmul(&Q_x, s_inv.as_le_bytes(), s_inv.nbits());

    let PQ_x = params.starting_curve.projective_difference(&P_x, &Q_x);

    BasisX::from_points(&P_x, &Q_x, &PQ_x)
}

#[rstest]
fn method1_produces_points_on_curve(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis = method1(&params, (&s, &s_inv));

    assert_eq!(
        params.starting_curve.is_on_curve(&basis.P.x()),
        SUCCESS_RETVAL,
        "Method 1 results in an x(P) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.Q.x()),
        SUCCESS_RETVAL,
        "Method 1 results in an x(Q) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.PQ.x()),
        SUCCESS_RETVAL,
        "Method 1 results in an x(P-Q) that is not on the curve",
    );
}

#[rstest]
fn method2_produces_points_on_curve(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis = method2(&params, (&s, &s_inv));

    assert_eq!(
        params.starting_curve.is_on_curve(&basis.P.x()),
        SUCCESS_RETVAL,
        "Method 2 results in an x(P) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.Q.x()),
        SUCCESS_RETVAL,
        "Method 2 results in an x(Q) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.PQ.x()),
        SUCCESS_RETVAL,
        "Method 2 results in an x(P-Q) that is not on the curve",
    );
}

#[rstest]
fn method3_produces_points_on_curve(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis = method3(&params, (&s, &s_inv));

    assert_eq!(
        params.starting_curve.is_on_curve(&basis.P.x()),
        SUCCESS_RETVAL,
        "Method 3 results in an x(P) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.Q.x()),
        SUCCESS_RETVAL,
        "Method 3 results in an x(Q) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.PQ.x()),
        SUCCESS_RETVAL,
        "Method 3 results in an x(P-Q) that is not on the curve",
    );
}

#[rstest]
fn projective_difference_method_produces_points_on_curve(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis = projective_difference_method(&params, (&s, &s_inv));

    assert_eq!(
        params.starting_curve.is_on_curve(&basis.P.x()),
        SUCCESS_RETVAL,
        "Projective Difference Method results in an x(P) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.Q.x()),
        SUCCESS_RETVAL,
        "Projective Difference Method results in an x(Q) that is not on the curve",
    );
    assert_eq!(
        params.starting_curve.is_on_curve(&basis.PQ.x()),
        SUCCESS_RETVAL,
        "Projective Difference Method results in an x(P-Q) that is not on the curve",
    );
}

#[rstest]
#[ignore]
fn method1_equals_method2(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis_method1 = method1(&params, (&s, &s_inv));
    let basis_method2 = method2(&params, (&s, &s_inv));

    assert_eq!(
        basis_method1.P.equals(&basis_method2.P),
        SUCCESS_RETVAL,
        "Method 1's\nx(P) = {}\ndoes not equal Method 2's\nx(P) = {}",
        basis_method1.P,
        basis_method2.P,
    );
    assert_eq!(
        basis_method1.Q.equals(&basis_method2.Q),
        SUCCESS_RETVAL,
        "Method 1's\nx(Q) = {}\ndoes not equal Method 2's\nx(Q) = {}",
        basis_method1.Q,
        basis_method2.Q,
    );
    assert_eq!(
        basis_method1.PQ.equals(&basis_method2.PQ),
        SUCCESS_RETVAL,
        "Method 1's\nx(P-Q) = {}\ndoes not equal Method 2's\nx(P-Q) = {}",
        basis_method1.PQ,
        basis_method2.PQ,
    );
}

#[rstest]
#[ignore]
fn method1_equals_method3(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis_method1 = method1(&params, (&s, &s_inv));
    let basis_method3 = method3(&params, (&s, &s_inv));

    assert_eq!(
        basis_method1.P.equals(&basis_method3.P),
        SUCCESS_RETVAL,
        "Method 1's\nx(P) = {}\ndoes not equal Method 3's\nx(P) = {}",
        basis_method1.P,
        basis_method3.P,
    );
    assert_eq!(
        basis_method1.Q.equals(&basis_method3.Q),
        SUCCESS_RETVAL,
        "Method 1's\nx(Q) = {}\ndoes not equal Method 3's\nx(Q) = {}",
        basis_method1.Q,
        basis_method3.Q,
    );
    assert_eq!(
        basis_method1.PQ.equals(&basis_method3.PQ),
        SUCCESS_RETVAL,
        "Method 1's\nx(P-Q) = {}\ndoes not equal Method 3's\nx(P-Q) = {}",
        basis_method1.PQ,
        basis_method3.PQ,
    );
}

#[rstest]
#[ignore]
fn method1_equals_projective_diff(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis_method1 = method1(&params, (&s, &s_inv));
    let basis_projective_difference_method = projective_difference_method(&params, (&s, &s_inv));

    assert_eq!(
        basis_method1
            .P
            .equals(&basis_projective_difference_method.P),
        SUCCESS_RETVAL,
        "Method 1's\nx(P) = {}\ndoes not equal Projective Difference Method's\nx(P) = {}",
        basis_method1.P,
        basis_projective_difference_method.P,
    );
    assert_eq!(
        basis_method1
            .Q
            .equals(&basis_projective_difference_method.Q),
        SUCCESS_RETVAL,
        "Method 1's\nx(Q) = {}\ndoes not equal Projective Difference Method's\nx(Q) = {}",
        basis_method1.Q,
        basis_projective_difference_method.Q,
    );
    assert_eq!(
        basis_method1
            .PQ
            .equals(&basis_projective_difference_method.PQ),
        SUCCESS_RETVAL,
        "Method 1's\nx(P-Q) = {}\ndoes not equal Projective Difference Method's\nx(P-Q) = {}",
        basis_method1.PQ,
        basis_projective_difference_method.PQ,
    );
}

#[rstest]
#[ignore]
fn method2_equals_method3(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis_method2 = method2(&params, (&s, &s_inv));
    let basis_method3 = method3(&params, (&s, &s_inv));

    assert_eq!(
        basis_method2.P.equals(&basis_method3.P),
        SUCCESS_RETVAL,
        "Method 2's\nx(P) = {}\ndoes not equal Method 3's\nx(P) = {}",
        basis_method2.P,
        basis_method3.P,
    );
    assert_eq!(
        basis_method2.Q.equals(&basis_method3.Q),
        SUCCESS_RETVAL,
        "Method 2's\nx(Q) = {}\ndoes not equal Method 3's\nx(Q) = {}",
        basis_method2.Q,
        basis_method3.Q,
    );
    assert_eq!(
        basis_method2.PQ.equals(&basis_method3.PQ),
        SUCCESS_RETVAL,
        "Method 2's\nx(P-Q) = {}\ndoes not equal Method 3's\nx(P-Q) = {}",
        basis_method2.PQ,
        basis_method3.PQ,
    );
}

#[rstest]
#[ignore]
fn method2_equals_projective_diff(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis_method2 = method2(&params, (&s, &s_inv));
    let basis_projective_difference_method = projective_difference_method(&params, (&s, &s_inv));

    assert_eq!(
        basis_method2
            .P
            .equals(&basis_projective_difference_method.P),
        SUCCESS_RETVAL,
        "Method 2's\nx(P) = {}\ndoes not equal Projective Difference Method's\nx(P) = {}",
        basis_method2.P,
        basis_projective_difference_method.P,
    );
    assert_eq!(
        basis_method2
            .Q
            .equals(&basis_projective_difference_method.Q),
        SUCCESS_RETVAL,
        "Method 2's\nx(Q) = {}\ndoes not equal Projective Difference Method's\nx(Q) = {}",
        basis_method2.Q,
        basis_projective_difference_method.Q,
    );
    assert_eq!(
        basis_method2
            .PQ
            .equals(&basis_projective_difference_method.PQ),
        SUCCESS_RETVAL,
        "Method 2's\nx(P-Q) = {}\ndoes not equal Projective Difference Method's\nx(P-Q) = {}",
        basis_method2.PQ,
        basis_projective_difference_method.PQ,
    );
}

#[rstest]
#[ignore]
fn method3_equals_projective_diff(
    params: PublicParams<PokeFieldI>,
    #[from(scalars)] (s, s_inv): (BigNum, BigNum),
) {
    let basis_method3 = method3(&params, (&s, &s_inv));
    let basis_projective_difference_method = projective_difference_method(&params, (&s, &s_inv));

    assert_eq!(
        basis_method3
            .P
            .equals(&basis_projective_difference_method.P),
        SUCCESS_RETVAL,
        "Method 3's\nx(P) = {}\ndoes not equal Projective Difference Method's\nx(P) = {}",
        basis_method3.P,
        basis_projective_difference_method.P,
    );
    assert_eq!(
        basis_method3
            .Q
            .equals(&basis_projective_difference_method.Q),
        SUCCESS_RETVAL,
        "Method 3's\nx(Q) = {}\ndoes not equal Projective Difference Method's\nx(Q) = {}",
        basis_method3.Q,
        basis_projective_difference_method.Q,
    );
    assert_eq!(
        basis_method3
            .PQ
            .equals(&basis_projective_difference_method.PQ),
        SUCCESS_RETVAL,
        "Method 3's\nx(P-Q) = {}\ndoes not equal Projective Difference Method's\nx(P-Q) = {}",
        basis_method3.PQ,
        basis_projective_difference_method.PQ,
    );
}
