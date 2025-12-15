#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use isogeny::elliptic::basis::BasisX;
use num_bigint::BigUint;
use poke::{bn::BigNum, params::poke_i, rand::sample_random_unit_mod_prime_power};

fn scalar_multiplication_then_basis_reconstruction(c: &mut Criterion) {
    let params = poke_i::get_params();

    // Generate scalars by which to multiply basis points
    let s1 = sample_random_unit_mod_prime_power(2, &params.full_two_torsion_order);
    let s2 = sample_random_unit_mod_prime_power(2, &params.full_two_torsion_order);

    // Benchmark the different methods to reconstruct an x-only basis after multiplying the 2 points in it
    let mut group = c.benchmark_group("Multiply then reconstruct basis/POKÉ level I");
    group.bench_function(
        "Method 1: Lift x-only basis to P, Q -> Multiply P, Q by s, t -> Subtract [s]*P - [t]*Q -> Convert to PointX -> Create basis",
        |b| {
            b.iter(|| {
                let (P, Q) = params.starting_curve.lift_basis(&params.two_torsion_basis);

                let P = params.starting_curve.mul(&P, s1.as_le_bytes(), s1.nbits());
                let Q = params.starting_curve.mul(&Q, s2.as_le_bytes(), s2.nbits());

                let PQ = params.starting_curve.sub(&P, &Q);

                BasisX::from_points(&P.to_pointx(), &Q.to_pointx(), &PQ.to_pointx())
            })
        },
    );
    group.bench_function(
        "Method 2: x-only multiply P, Q by s, t -> Lift x([s]*P), x([t]*Q) to full point -> Subtract [s]*P - [t]*Q -> Convert to PointX -> Create basis",
        |b| b.iter(|| {
            let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();

            let P_x = params.starting_curve.xmul(&P_x, s1.as_le_bytes(), s1.nbits());
            let Q_x = params.starting_curve.xmul(&Q_x, s2.as_le_bytes(), s2.nbits());

            let (P, _) = params.starting_curve.lift_pointx(&P_x);
            let (Q, _) = params.starting_curve.lift_pointx(&Q_x);

            let PQ = params.starting_curve.sub(&P, &Q);

            BasisX::from_points(&P_x, &Q_x, &PQ.to_pointx())
        }),
    );
    group.bench_function(
        "Method 3: x-only multiply P, Q by s, t -> Compute x([s]*P - [t]*Q) using biscalar ladder -> Create basis",
        |b| {
            b.iter(|| {
                let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();

                let P_x = params.starting_curve.xmul(&P_x, s1.as_le_bytes(), s1.nbits());
                let Q_x = params.starting_curve.xmul(&Q_x, s2.as_le_bytes(), s2.nbits());

                let minus_s2 = BigUint::from_bytes_le(params.full_two_torsion_order.as_le_bytes()) - BigUint::from_bytes_le(s2.as_le_bytes());
                let minus_s2 = BigNum::new(&minus_s2.to_u64_digits());

                let PQ_x = params.starting_curve.ladder_biscalar(
                    &params.two_torsion_basis,
                    s1.as_le_bytes(),
                    minus_s2.as_le_bytes(),
                    s1.nbits(),
                    minus_s2.nbits(),
                );

                BasisX::from_points(&P_x, &Q_x, &PQ_x)
            })
        },
    );
}

criterion_group!(benches, scalar_multiplication_then_basis_reconstruction);
criterion_main!(benches);
