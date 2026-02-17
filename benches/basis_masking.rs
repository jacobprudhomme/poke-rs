#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point};
use poke::{
    bn::BigNum,
    params::{poke_i, poke_iii, poke_v},
    rand::{
        sample_random_invertible_matrix_mod_prime_power, sample_random_torsion_basis,
        sample_random_unit_mod_prime_power,
    },
};

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

// fn multiply_xonly_basis_by_same_scalar_projective_difference<
//     Fp2: Fp2Trait,
//     const NUM_WORDS: usize,
// >(
//     curve: &Curve<Fp2>,
//     basis: &BasisX<Fp2>,
//     s: &BigNum<NUM_WORDS>,
// ) -> BasisX<Fp2> {
//     let masked_xP = curve.xmul(&basis.P, &s.to_le_bytes(), s.nbits());
//     let masked_xPQ = curve.xmul(&basis.PQ, &s.to_le_bytes(), s.nbits());
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

fn mask_basis_by_same_scalar(c: &mut Criterion) {
    // Use largest parameters to see difference most clearly
    let params_i = poke_i::get_params();
    let params_iii = poke_iii::get_params();
    let params_v = poke_v::get_params();

    // Generate scalars by which to multiply basis points
    let s_i = sample_random_unit_mod_prime_power(
        params_i.three_torsion.base,
        &params_i.three_torsion.order,
    );
    let s_iii = sample_random_unit_mod_prime_power(
        params_iii.three_torsion.base,
        &params_iii.three_torsion.order,
    );
    let s_v = sample_random_unit_mod_prime_power(
        params_v.three_torsion.base,
        &params_v.three_torsion.order,
    );

    // Generate random bases of given order
    let (P_i, Q_i, _) = sample_random_torsion_basis(
        &params_i.starting_curve,
        &[params_i.three_torsion.base],
        &params_i.three_torsion.order,
        &params_i.cofactor,
    );
    let xP_i = P_i.to_pointx();
    let xQ_i = Q_i.to_pointx();
    let xPQ_i = params_i.starting_curve.sub(&P_i, &Q_i).to_pointx();

    let (P_iii, Q_iii, _) = sample_random_torsion_basis(
        &params_iii.starting_curve,
        &[params_iii.three_torsion.base],
        &params_iii.three_torsion.order,
        &params_iii.cofactor,
    );
    let xP_iii = P_iii.to_pointx();
    let xQ_iii = Q_iii.to_pointx();
    let xPQ_iii = params_iii.starting_curve.sub(&P_iii, &Q_iii).to_pointx();

    let (P_v, Q_v, _) = sample_random_torsion_basis(
        &params_v.starting_curve,
        &[params_v.three_torsion.base],
        &params_v.three_torsion.order,
        &params_v.cofactor,
    );
    let xP_v = P_v.to_pointx();
    let xQ_v = Q_v.to_pointx();
    let xPQ_v = params_v.starting_curve.sub(&P_v, &Q_v).to_pointx();

    let mut lifted_basis_group =
        c.benchmark_group("Method 1: Full point multiplication on lifted basis");
    lifted_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            let basis = params_i.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalars(&params_i.starting_curve, &basis, &s_i, &s_i)
        })
    });
    lifted_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            let basis = params_iii.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalars(&params_iii.starting_curve, &basis, &s_iii, &s_iii)
        })
    });
    lifted_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            let basis = params_v.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalars(&params_v.starting_curve, &basis, &s_v, &s_v)
        })
    });
    lifted_basis_group.finish();

    let mut lifted_normalized_basis_group =
        c.benchmark_group("Method 2: Full point multiplication on lifted, normalized basis");
    lifted_normalized_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis =
                params_i
                    .starting_curve
                    .lift_basis_normalised(&xP_i.x(), &xQ_i.x(), &xPQ_i.x());
            multiply_basis_by_scalars(&params_i.starting_curve, &basis, &s_i, &s_i)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = params_iii.starting_curve.lift_basis_normalised(
                &xP_iii.x(),
                &xQ_iii.x(),
                &xPQ_iii.x(),
            );
            multiply_basis_by_scalars(&params_iii.starting_curve, &basis, &s_iii, &s_iii)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis =
                params_v
                    .starting_curve
                    .lift_basis_normalised(&xP_v.x(), &xQ_v.x(), &xPQ_v.x());
            multiply_basis_by_scalars(&params_v.starting_curve, &basis, &s_v, &s_v)
        })
    });
    lifted_normalized_basis_group.finish();

    let mut xmul_group = c.benchmark_group("Method 3: x-only point multiplication");
    xmul_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            multiply_xonly_basis_by_same_scalar_xmul(&params_i.starting_curve, &basis, &s_i)
        })
    });
    xmul_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            multiply_xonly_basis_by_same_scalar_xmul(&params_iii.starting_curve, &basis, &s_iii)
        })
    });
    xmul_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            multiply_xonly_basis_by_same_scalar_xmul(&params_v.starting_curve, &basis, &s_v)
        })
    });
    xmul_group.finish();

    // let mut projective_difference_group =
    //     c.benchmark_group("Method 4: x-only point multiplication and projective difference");
    // projective_difference_group.bench_function("POKÉ level I", |b| {
    //     b.iter(|| {
    //         let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
    //         multiply_xonly_basis_by_same_scalar_projective_difference(
    //             &params_i.starting_curve,
    //             &basis,
    //             &s_i,
    //         )
    //     })
    // });
    // projective_difference_group.bench_function("POKÉ level III", |b| {
    //     b.iter(|| {
    //         let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
    //         multiply_xonly_basis_by_same_scalar_projective_difference(
    //             &params_iii.starting_curve,
    //             &basis,
    //             &s_iii,
    //         )
    //     })
    // });
    // projective_difference_group.bench_function("POKÉ level V", |b| {
    //     b.iter(|| {
    //         let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
    //         multiply_xonly_basis_by_same_scalar_projective_difference(
    //             &params_v.starting_curve,
    //             &basis,
    //             &s_v,
    //         )
    //     })
    // });
    // projective_difference_group.finish();
}

fn mask_basis_by_different_scalars(c: &mut Criterion) {
    // Use largest parameters to see difference most clearly
    let params_i = poke_i::get_params();
    let params_iii = poke_iii::get_params();
    let params_v = poke_v::get_params();

    // Generate scalars by which to multiply basis points
    let s1_i = sample_random_unit_mod_prime_power(
        params_i.three_torsion.base,
        &params_i.three_torsion.order,
    );
    let s2_i = sample_random_unit_mod_prime_power(
        params_i.three_torsion.base,
        &params_i.three_torsion.order,
    );
    let s1_iii = sample_random_unit_mod_prime_power(
        params_iii.three_torsion.base,
        &params_iii.three_torsion.order,
    );
    let s2_iii = sample_random_unit_mod_prime_power(
        params_iii.three_torsion.base,
        &params_iii.three_torsion.order,
    );
    let s1_v = sample_random_unit_mod_prime_power(
        params_v.three_torsion.base,
        &params_v.three_torsion.order,
    );
    let s2_v = sample_random_unit_mod_prime_power(
        params_v.three_torsion.base,
        &params_v.three_torsion.order,
    );

    // Generate random bases of given order
    let (P_i, Q_i, _) = sample_random_torsion_basis(
        &params_i.starting_curve,
        &[params_i.three_torsion.base],
        &params_i.three_torsion.order,
        &params_i.cofactor,
    );
    let xP_i = P_i.to_pointx();
    let xQ_i = Q_i.to_pointx();
    let xPQ_i = params_i.starting_curve.sub(&P_i, &Q_i).to_pointx();

    let (P_iii, Q_iii, _) = sample_random_torsion_basis(
        &params_iii.starting_curve,
        &[params_iii.three_torsion.base],
        &params_iii.three_torsion.order,
        &params_iii.cofactor,
    );
    let xP_iii = P_iii.to_pointx();
    let xQ_iii = Q_iii.to_pointx();
    let xPQ_iii = params_iii.starting_curve.sub(&P_iii, &Q_iii).to_pointx();

    let (P_v, Q_v, _) = sample_random_torsion_basis(
        &params_v.starting_curve,
        &[params_v.three_torsion.base],
        &params_v.three_torsion.order,
        &params_v.cofactor,
    );
    let xP_v = P_v.to_pointx();
    let xQ_v = Q_v.to_pointx();
    let xPQ_v = params_v.starting_curve.sub(&P_v, &Q_v).to_pointx();

    let mut lifted_basis_group =
        c.benchmark_group("Method 1: Full point multiplication on lifted basis");
    lifted_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            let basis = params_i.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalars(&params_i.starting_curve, &basis, &s1_i, &s2_i)
        })
    });
    lifted_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            let basis = params_iii.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalars(&params_iii.starting_curve, &basis, &s1_iii, &s2_iii)
        })
    });
    lifted_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            let basis = params_v.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalars(&params_v.starting_curve, &basis, &s1_v, &s2_v)
        })
    });
    lifted_basis_group.finish();

    let mut lifted_normalized_basis_group =
        c.benchmark_group("Method 2: Full point multiplication on lifted, normalized basis");
    lifted_normalized_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis =
                params_i
                    .starting_curve
                    .lift_basis_normalised(&xP_i.x(), &xQ_i.x(), &xPQ_i.x());
            multiply_basis_by_scalars(&params_i.starting_curve, &basis, &s1_i, &s2_i)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = params_iii.starting_curve.lift_basis_normalised(
                &xP_iii.x(),
                &xQ_iii.x(),
                &xPQ_iii.x(),
            );
            multiply_basis_by_scalars(&params_iii.starting_curve, &basis, &s1_iii, &s2_iii)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis =
                params_v
                    .starting_curve
                    .lift_basis_normalised(&xP_v.x(), &xQ_v.x(), &xPQ_v.x());
            multiply_basis_by_scalars(&params_v.starting_curve, &basis, &s1_v, &s2_v)
        })
    });
    lifted_normalized_basis_group.finish();

    let mut flipped_Q_group =
        c.benchmark_group("Method 3: Biscalar ladder using an x-only basis of (P,-Q)");
    flipped_Q_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            multiply_xonly_basis_by_scalars_flipped_Q_basis(
                &params_i.starting_curve,
                &basis,
                &s1_i,
                &s2_i,
            )
        })
    });
    flipped_Q_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            multiply_xonly_basis_by_scalars_flipped_Q_basis(
                &params_iii.starting_curve,
                &basis,
                &s1_iii,
                &s2_iii,
            )
        })
    });
    flipped_Q_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            multiply_xonly_basis_by_scalars_flipped_Q_basis(
                &params_v.starting_curve,
                &basis,
                &s1_v,
                &s2_v,
            )
        })
    });
    flipped_Q_group.finish();

    let mut neg_s2_group = c.benchmark_group("Method 4: Biscalar ladder using [-s2]");
    neg_s2_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            multiply_xonly_basis_by_scalars_negate_second_scalar(
                &params_i.starting_curve,
                &basis,
                &s1_i,
                &s2_i,
                &params_i.three_torsion.order,
            )
        })
    });
    neg_s2_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            multiply_xonly_basis_by_scalars_negate_second_scalar(
                &params_iii.starting_curve,
                &basis,
                &s1_iii,
                &s2_iii,
                &params_iii.three_torsion.order,
            )
        })
    });
    neg_s2_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            multiply_xonly_basis_by_scalars_negate_second_scalar(
                &params_v.starting_curve,
                &basis,
                &s1_v,
                &s2_v,
                &params_v.three_torsion.order,
            )
        })
    });
    neg_s2_group.finish();

    let mut inv_s1_group = c.benchmark_group("Method 5: Biscalar ladder using [-s2/s1]");
    inv_s1_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            multiply_xonly_basis_by_scalars_using_invert_first_scalar(
                &params_i.starting_curve,
                &basis,
                &s1_i,
                &s2_i,
                &params_i.three_torsion.order,
            )
        })
    });
    inv_s1_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            multiply_xonly_basis_by_scalars_using_invert_first_scalar(
                &params_iii.starting_curve,
                &basis,
                &s1_iii,
                &s2_iii,
                &params_iii.three_torsion.order,
            )
        })
    });
    inv_s1_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            multiply_xonly_basis_by_scalars_using_invert_first_scalar(
                &params_v.starting_curve,
                &basis,
                &s1_v,
                &s2_v,
                &params_v.three_torsion.order,
            )
        })
    });
    inv_s1_group.finish();
}

fn mask_basis_by_scalar_matrix(c: &mut Criterion) {
    // Use largest parameters to see difference most clearly
    let params_i = poke_i::get_params();
    let params_iii = poke_iii::get_params();
    let params_v = poke_v::get_params();

    // Generate scalars by which to multiply basis points
    let S_i = sample_random_invertible_matrix_mod_prime_power(
        params_i.three_torsion.base,
        &params_i.three_torsion.order,
    );
    let S_iii = sample_random_invertible_matrix_mod_prime_power(
        params_iii.three_torsion.base,
        &params_iii.three_torsion.order,
    );
    let S_v = sample_random_invertible_matrix_mod_prime_power(
        params_v.three_torsion.base,
        &params_v.three_torsion.order,
    );

    // Generate random bases of given order
    let (P_i, Q_i, _) = sample_random_torsion_basis(
        &params_i.starting_curve,
        &[params_i.three_torsion.base],
        &params_i.three_torsion.order,
        &params_i.cofactor,
    );
    let xP_i = P_i.to_pointx();
    let xQ_i = Q_i.to_pointx();
    let xPQ_i = params_i.starting_curve.sub(&P_i, &Q_i).to_pointx();

    let (P_iii, Q_iii, _) = sample_random_torsion_basis(
        &params_iii.starting_curve,
        &[params_iii.three_torsion.base],
        &params_iii.three_torsion.order,
        &params_iii.cofactor,
    );
    let xP_iii = P_iii.to_pointx();
    let xQ_iii = Q_iii.to_pointx();
    let xPQ_iii = params_iii.starting_curve.sub(&P_iii, &Q_iii).to_pointx();

    let (P_v, Q_v, _) = sample_random_torsion_basis(
        &params_v.starting_curve,
        &[params_v.three_torsion.base],
        &params_v.three_torsion.order,
        &params_v.cofactor,
    );
    let xP_v = P_v.to_pointx();
    let xQ_v = Q_v.to_pointx();
    let xPQ_v = params_v.starting_curve.sub(&P_v, &Q_v).to_pointx();

    let mut lifted_basis_group =
        c.benchmark_group("Method 1: Full point multiplication on lifted basis");
    lifted_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            let basis = params_i.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalar_matrix(&params_i.starting_curve, &basis, &S_i)
        })
    });
    lifted_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            let basis = params_iii.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalar_matrix(&params_iii.starting_curve, &basis, &S_iii)
        })
    });
    lifted_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            let basis = params_v.starting_curve.lift_basis(&basisx);
            multiply_basis_by_scalar_matrix(&params_v.starting_curve, &basis, &S_v)
        })
    });
    lifted_basis_group.finish();

    let mut lifted_normalized_basis_group =
        c.benchmark_group("Method 2: Full point multiplication on lifted, normalized basis");
    lifted_normalized_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis =
                params_i
                    .starting_curve
                    .lift_basis_normalised(&xP_i.x(), &xQ_i.x(), &xPQ_i.x());
            multiply_basis_by_scalar_matrix(&params_i.starting_curve, &basis, &S_i)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = params_iii.starting_curve.lift_basis_normalised(
                &xP_iii.x(),
                &xQ_iii.x(),
                &xPQ_iii.x(),
            );
            multiply_basis_by_scalar_matrix(&params_iii.starting_curve, &basis, &S_iii)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis =
                params_v
                    .starting_curve
                    .lift_basis_normalised(&xP_v.x(), &xQ_v.x(), &xPQ_v.x());
            multiply_basis_by_scalar_matrix(&params_v.starting_curve, &basis, &S_v)
        })
    });
    lifted_normalized_basis_group.finish();

    let mut scalar_difference_group =
        c.benchmark_group("Method 3: Biscalar ladder using column-wise scalar differences");
    scalar_difference_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            multiply_xonly_basis_by_scalar_matrix_scalar_difference(
                &params_i.starting_curve,
                &basis,
                &S_i,
                &params_i.three_torsion.order,
            )
        })
    });
    scalar_difference_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            multiply_xonly_basis_by_scalar_matrix_scalar_difference(
                &params_iii.starting_curve,
                &basis,
                &S_iii,
                &params_iii.three_torsion.order,
            )
        })
    });
    scalar_difference_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            multiply_xonly_basis_by_scalar_matrix_scalar_difference(
                &params_v.starting_curve,
                &basis,
                &S_v,
                &params_v.three_torsion.order,
            )
        })
    });
    scalar_difference_group.finish();
}

fn special_case_mask_basis_by_scalar_matrix_and_keep_xP_xQ_only(c: &mut Criterion) {
    // Use largest parameters to see difference most clearly
    let params_i = poke_i::get_params();
    let params_iii = poke_iii::get_params();
    let params_v = poke_v::get_params();

    // Generate scalars by which to multiply basis points
    let S_i = sample_random_invertible_matrix_mod_prime_power(
        params_i.three_torsion.base,
        &params_i.three_torsion.order,
    );
    let S_iii = sample_random_invertible_matrix_mod_prime_power(
        params_iii.three_torsion.base,
        &params_iii.three_torsion.order,
    );
    let S_v = sample_random_invertible_matrix_mod_prime_power(
        params_v.three_torsion.base,
        &params_v.three_torsion.order,
    );

    // Generate random bases of given order
    let (P_i, Q_i, _) = sample_random_torsion_basis(
        &params_i.starting_curve,
        &[params_i.three_torsion.base],
        &params_i.three_torsion.order,
        &params_i.cofactor,
    );
    let xP_i = P_i.to_pointx();
    let xQ_i = Q_i.to_pointx();
    let xPQ_i = params_i.starting_curve.sub(&P_i, &Q_i).to_pointx();

    let (P_iii, Q_iii, _) = sample_random_torsion_basis(
        &params_iii.starting_curve,
        &[params_iii.three_torsion.base],
        &params_iii.three_torsion.order,
        &params_iii.cofactor,
    );
    let xP_iii = P_iii.to_pointx();
    let xQ_iii = Q_iii.to_pointx();
    let xPQ_iii = params_iii.starting_curve.sub(&P_iii, &Q_iii).to_pointx();

    let (P_v, Q_v, _) = sample_random_torsion_basis(
        &params_v.starting_curve,
        &[params_v.three_torsion.base],
        &params_v.three_torsion.order,
        &params_v.cofactor,
    );
    let xP_v = P_v.to_pointx();
    let xQ_v = Q_v.to_pointx();
    let xPQ_v = params_v.starting_curve.sub(&P_v, &Q_v).to_pointx();

    let mut lifted_basis_group =
        c.benchmark_group("Method 1: Full point multiplication on lifted basis");
    lifted_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            let basis = params_i.starting_curve.lift_basis(&basisx);
            special_case_multiply_basis_by_scalar_matrix(&params_i.starting_curve, &basis, &S_i)
        })
    });
    lifted_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            let basis = params_iii.starting_curve.lift_basis(&basisx);
            special_case_multiply_basis_by_scalar_matrix(&params_iii.starting_curve, &basis, &S_iii)
        })
    });
    lifted_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basisx = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            let basis = params_v.starting_curve.lift_basis(&basisx);
            special_case_multiply_basis_by_scalar_matrix(&params_v.starting_curve, &basis, &S_v)
        })
    });
    lifted_basis_group.finish();

    let mut lifted_normalized_basis_group =
        c.benchmark_group("Method 2: Full point multiplication on lifted, normalized basis");
    lifted_normalized_basis_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis =
                params_i
                    .starting_curve
                    .lift_basis_normalised(&xP_i.x(), &xQ_i.x(), &xPQ_i.x());
            special_case_multiply_basis_by_scalar_matrix(&params_i.starting_curve, &basis, &S_i)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = params_iii.starting_curve.lift_basis_normalised(
                &xP_iii.x(),
                &xQ_iii.x(),
                &xPQ_iii.x(),
            );
            special_case_multiply_basis_by_scalar_matrix(&params_iii.starting_curve, &basis, &S_iii)
        })
    });
    lifted_normalized_basis_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis =
                params_v
                    .starting_curve
                    .lift_basis_normalised(&xP_v.x(), &xQ_v.x(), &xPQ_v.x());
            special_case_multiply_basis_by_scalar_matrix(&params_v.starting_curve, &basis, &S_v)
        })
    });
    lifted_normalized_basis_group.finish();

    let mut biscalar_ladder_group = c.benchmark_group("Method 3: Biscalar ladder");
    biscalar_ladder_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_i, &xQ_i, &xPQ_i);
            special_case_multiply_xonly_basis_by_scalar_matrix(
                &params_i.starting_curve,
                &basis,
                &S_i,
            )
        })
    });
    biscalar_ladder_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_iii, &xQ_iii, &xPQ_iii);
            special_case_multiply_xonly_basis_by_scalar_matrix(
                &params_iii.starting_curve,
                &basis,
                &S_iii,
            )
        })
    });
    biscalar_ladder_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            let basis = BasisX::from_points(&xP_v, &xQ_v, &xPQ_v);
            special_case_multiply_xonly_basis_by_scalar_matrix(
                &params_v.starting_curve,
                &basis,
                &S_v,
            )
        })
    });
    biscalar_ladder_group.finish();
}

criterion_group!(
    benches,
    mask_basis_by_same_scalar,
    mask_basis_by_different_scalars,
    mask_basis_by_scalar_matrix,
    special_case_mask_basis_by_scalar_matrix_and_keep_xP_xQ_only,
);
criterion_main!(benches);
