#![allow(non_snake_case)]

use criterion::{criterion_group, criterion_main, Criterion};
use isogeny::elliptic::basis::BasisX;
use num_bigint::{BigUint, RandBigInt as _};
use poke::poke_i::create_poke_i_params;

fn basis_reconstruction(c: &mut Criterion) {
    let params = create_poke_i_params();

    // Generate scalars by which to multiply basis points
    let mut rng = rand::thread_rng();
    let ONE = BigUint::from(1u8);
    let Z_two_torsion_order = BigUint::from(2u8).pow(
        params
            .two_torsion_exp
            .try_into()
            .expect("Exponent of the 2-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let mut s = rng.gen_biguint_range(&ONE, &Z_two_torsion_order);
    let mut s_inv = s.modinv(&Z_two_torsion_order);
    while let None = s_inv {
        s = rng.gen_biguint_range(&ONE, &Z_two_torsion_order);
        s_inv = s.modinv(&Z_two_torsion_order);
    }
    let Some(s_inv) = s_inv else {
        unreachable!();
    };
    let s_bitsize =
        s.bits().try_into().expect("Size in bits of the scalar s is too big to fit in a usize (we do not ever expect this to happen)");
    let s_inv_bitsize =
        s_inv.bits().try_into().expect("Size in bits of the scalar 1/s is too big to fit in a usize (we do not ever expect this to happen)");
    let s = s.to_bytes_le();
    let s_inv = s_inv.to_bytes_le();

    // Multiply points by these scalars
    let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();
    let P_x = params.starting_curve.xmul(&P_x, &s, s_bitsize);
    let Q_x = params.starting_curve.xmul(&Q_x, &s_inv, s_inv_bitsize);

    // Precompute values and compare to make sure the different methods give the same answer, before proceeding
    let (P, _) = params.starting_curve.lift_pointx(&P_x);
    let (Q, _) = params.starting_curve.lift_pointx(&Q_x);
    let PQ_method1 = params.starting_curve.sub(&P, &Q);

    let PQ_x = params.starting_curve.projective_difference(&P_x, &Q_x);
    let (PQ_method2, _) = params.starting_curve.lift_pointx(&PQ_x);

    assert_eq!(PQ_method1.equals(&PQ_method2), 0xffffffff);

    let mut group_poke_i = c.benchmark_group("POKÉ level I");
    group_poke_i.bench_function(
        "Lift x(P), x(Q) to full point -> Subtract P - Q -> Convert to PointX -> Create basis",
        |b| b.iter(|| {
            let (P, _) = params.starting_curve.lift_pointx(&P_x);
            let (Q, _) = params.starting_curve.lift_pointx(&Q_x);
            let PQ = params.starting_curve.sub(&P, &Q);
            BasisX::from_points(&P_x, &Q_x, &PQ.to_pointx())
        }),
    );
    group_poke_i.bench_function(
        "Projective subtraction between P and Q -> Create basis",
        |b| b.iter(|| {
            let PQ_x = params.starting_curve.projective_difference(&P_x, &Q_x);
            BasisX::from_points(&P_x, &Q_x, &PQ_x)
        }),
    );
}

criterion_group!(benches, basis_reconstruction);
criterion_main!(benches);
