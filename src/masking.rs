use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point};

use crate::bn::BigNum;

pub fn mask_basis_by_same_scalar<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    (P, Q): &(Point<Fp2>, Point<Fp2>),
    scalar: &BigNum<NUM_WORDS>,
) -> (Point<Fp2>, Point<Fp2>) {
    let scalar_bitlen = scalar.nbits();
    let scalar = scalar.to_le_bytes();

    let masked_P = curve.mul(P, &scalar, scalar_bitlen);
    let masked_Q = curve.mul(Q, &scalar, scalar_bitlen);

    (masked_P, masked_Q)
}

pub fn mask_basisx_by_same_scalar<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    scalar: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let (masked_P, masked_Q) = mask_basis_by_same_scalar(curve, &curve.lift_basis(basis), scalar);
    let masked_PQ = curve.sub(&masked_P, &masked_Q);

    BasisX::from_points(
        &masked_P.to_pointx(),
        &masked_Q.to_pointx(),
        &masked_PQ.to_pointx(),
    )
}

pub fn mask_basisx_by_diagonal_scalars<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    scalar1: &BigNum<NUM_WORDS>,
    scalar2: &BigNum<NUM_WORDS>,
) -> BasisX<Fp2> {
    let (P, Q) = curve.lift_basis(basis);

    let masked_P = curve.mul(&P, &scalar1.to_le_bytes(), scalar1.nbits());
    let masked_Q = curve.mul(&Q, &scalar2.to_le_bytes(), scalar2.nbits());
    let masked_PQ = curve.sub(&masked_P, &masked_Q);

    BasisX::from_points(
        &masked_P.to_pointx(),
        &masked_Q.to_pointx(),
        &masked_PQ.to_pointx(),
    )
}

pub fn mask_basisx_by_diagonal_scalars_points_only<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    scalar1: &BigNum<NUM_WORDS>,
    scalar2: &BigNum<NUM_WORDS>,
) -> (Point<Fp2>, Point<Fp2>) {
    let (P, Q) = curve.lift_basis(basis);

    let masked_P = curve.mul(&P, &scalar1.to_le_bytes(), scalar1.nbits());
    let masked_Q = curve.mul(&Q, &scalar2.to_le_bytes(), scalar2.nbits());

    (masked_P, masked_Q)
}

pub fn mask_basisx_by_scalar_matrix<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    matrix: &[[BigNum<NUM_WORDS>; 2]; 2],
) -> BasisX<Fp2> {
    let (P, Q) = curve.lift_basis(basis);

    let masked_P = curve.add(
        &curve.mul(&P, &matrix[0][0].to_le_bytes(), matrix[0][0].nbits()),
        &curve.mul(&Q, &matrix[0][1].to_le_bytes(), matrix[0][1].nbits()),
    );
    let masked_Q = curve.add(
        &curve.mul(&P, &matrix[1][0].to_le_bytes(), matrix[1][0].nbits()),
        &curve.mul(&Q, &matrix[1][1].to_le_bytes(), matrix[1][1].nbits()),
    );
    let masked_PQ = curve.sub(&masked_P, &masked_Q);

    BasisX::from_points(
        &masked_P.to_pointx(),
        &masked_Q.to_pointx(),
        &masked_PQ.to_pointx(),
    )
}

pub fn mask_basisx_by_scalar_matrix_pointx_only<Fp2: Fp2Trait, const NUM_WORDS: usize>(
    curve: &Curve<Fp2>,
    basis: &BasisX<Fp2>,
    matrix: &[[BigNum<NUM_WORDS>; 2]; 2],
) -> (PointX<Fp2>, PointX<Fp2>) {
    let masked_xP = curve.ladder_biscalar(
        basis,
        &matrix[0][0].to_le_bytes(),
        &matrix[0][1].to_le_bytes(),
        matrix[0][0].nbits(),
        matrix[0][1].nbits(),
    );
    let masked_xQ = curve.ladder_biscalar(
        basis,
        &matrix[1][0].to_le_bytes(),
        &matrix[1][1].to_le_bytes(),
        matrix[1][0].nbits(),
        matrix[1][1].nbits(),
    );

    (masked_xP, masked_xQ)
}
