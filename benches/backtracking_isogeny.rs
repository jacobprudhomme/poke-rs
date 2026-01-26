#![allow(non_snake_case)]

use core::mem;

use criterion::{Criterion, criterion_group, criterion_main};
use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point},
    polynomial_ring::poly::Polynomial,
};
use poke::params::{poke_i, poke_iii, poke_v};

fn apply_backtracking_isogeny_three_chain_then_prime_power<Fp2: Fp2Trait>(
    curve: &Curve<Fp2>,
    K3: &Point<Fp2>,
    three_isogeny_length: usize,
    K5: &Point<Fp2>,
    five_isogeny_length: usize,
    basis: &BasisX<Fp2>,
) -> (Curve<Fp2>, BasisX<Fp2>) {
    let mut basis_and_K5 = [PointX::INFINITY; 4];
    basis_and_K5[..3].copy_from_slice(&basis.to_array());
    basis_and_K5[3] = K5.to_pointx();
    let (codomain, _) =
        curve.three_isogeny_chain(&K3.to_pointx(), three_isogeny_length, &mut basis_and_K5);

    let K5_img_under_backtracking_isogeny = mem::replace(&mut basis_and_K5[3], PointX::INFINITY);
    let mut basis = &mut basis_and_K5[..3];
    let codomain = codomain.velu_prime_power_isogeny::<Polynomial<Fp2>>(
        &K5_img_under_backtracking_isogeny,
        5,
        five_isogeny_length,
        &mut basis,
    );

    (codomain, BasisX::from_slice(&basis))
}

fn apply_backtracking_isogeny_composite<Fp2: Fp2Trait>(
    curve: &Curve<Fp2>,
    K3: &Point<Fp2>,
    three_isogeny_length: usize,
    K5: &Point<Fp2>,
    five_isogeny_length: usize,
    basis: &BasisX<Fp2>,
) -> (Curve<Fp2>, BasisX<Fp2>) {
    let mut basis = basis.to_array();
    let K35 = curve.add(&K3, &K5);
    let codomain = curve.velu_composite_isogeny::<Polynomial<Fp2>>(
        &K35.to_pointx(),
        &[(3, three_isogeny_length), (5, five_isogeny_length)],
        &mut basis,
    );

    (codomain, BasisX::from_slice(&basis))
}

fn backtracking_isogeny(c: &mut Criterion) {
    let params_i = poke_i::get_params();
    let params_iii = poke_iii::get_params();
    let params_v = poke_v::get_params();

    // Lift a single basis x-coordinate to a full point (doesn't matter the sign or which point, it will always have full order)
    let (K3_i, _) = params_i
        .starting_curve
        .lift_pointx(&params_i.three_torsion_basis.P);
    let (K5_i, _) = params_i
        .starting_curve
        .lift_pointx(&params_i.five_torsion_basis.P);
    let (K3_iii, _) = params_iii
        .starting_curve
        .lift_pointx(&params_iii.three_torsion_basis.P);
    let (K5_iii, _) = params_iii
        .starting_curve
        .lift_pointx(&params_iii.five_torsion_basis.P);
    let (K3_v, _) = params_v
        .starting_curve
        .lift_pointx(&params_v.three_torsion_basis.P);
    let (K5_v, _) = params_v
        .starting_curve
        .lift_pointx(&params_v.five_torsion_basis.P);

    let mut three_then_prime_power_group = c.benchmark_group(
        "Method 1: Applying a specialized 3-isogeny chain, then a degree-5^c prime-power isogeny",
    );
    three_then_prime_power_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            apply_backtracking_isogeny_three_chain_then_prime_power(
                &params_i.starting_curve,
                &K3_i,
                params_i.three_torsion.exp,
                &K5_i,
                params_i.five_torsion.exp,
                &params_i.two_torsion_basis,
            )
        })
    });
    three_then_prime_power_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            apply_backtracking_isogeny_three_chain_then_prime_power(
                &params_iii.starting_curve,
                &K3_iii,
                params_iii.three_torsion.exp,
                &K5_iii,
                params_iii.five_torsion.exp,
                &params_iii.two_torsion_basis,
            )
        })
    });
    three_then_prime_power_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            apply_backtracking_isogeny_three_chain_then_prime_power(
                &params_v.starting_curve,
                &K3_v,
                params_v.three_torsion.exp,
                &K5_v,
                params_v.five_torsion.exp,
                &params_v.two_torsion_basis,
            )
        })
    });
    three_then_prime_power_group.finish();

    let mut composite_group =
        c.benchmark_group("Method 2: Applying a composite, degree-(3^b*5^c) isogeny all at once");
    composite_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            apply_backtracking_isogeny_composite(
                &params_i.starting_curve,
                &K3_i,
                params_i.three_torsion.exp,
                &K5_i,
                params_i.five_torsion.exp,
                &params_i.two_torsion_basis,
            )
        })
    });
    composite_group.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            apply_backtracking_isogeny_composite(
                &params_iii.starting_curve,
                &K3_iii,
                params_iii.three_torsion.exp,
                &K5_iii,
                params_iii.five_torsion.exp,
                &params_iii.two_torsion_basis,
            )
        })
    });
    composite_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            apply_backtracking_isogeny_composite(
                &params_v.starting_curve,
                &K3_v,
                params_v.three_torsion.exp,
                &K5_v,
                params_v.five_torsion.exp,
                &params_v.two_torsion_basis,
            )
        })
    });
    composite_group.finish();
}

criterion_group!(benches, backtracking_isogeny);
criterion_main!(benches);
