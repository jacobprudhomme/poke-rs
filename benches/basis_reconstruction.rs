#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use isogeny::elliptic::basis::BasisX;
use num_bigint::{BigUint, RandBigInt as _};
use poke::{fields::PokeFieldIBase, poke_i::create_poke_i_params};

fn scalar_multiplication_then_basis_reconstruction(c: &mut Criterion) {
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

    // Benchmark the different methods to reconstruct an x-only basis after multiplying the 2 points in it
    let mut group = c.benchmark_group("Multiply then reconstruct basis/POKÉ level I");
    group.bench_function(
        "Method 1: Lift x-only basis to P, Q -> Multiply P, Q by s, t -> Subtract [s]*P - [t]*Q -> Convert to PointX -> Create basis",
        |b| {
            b.iter(|| {
                let (P, Q) = params.starting_curve.lift_basis(&params.two_torsion_basis);

                let P = params.starting_curve.mul(&P, &s, s_bitsize);
                let Q = params.starting_curve.mul(&Q, &s_inv, s_inv_bitsize);

                let PQ = params.starting_curve.sub(&P, &Q);

                BasisX::from_points(&P.to_pointx(), &Q.to_pointx(), &PQ.to_pointx())
            })
        },
    );
    group.bench_function(
        "Method 2: x-only multiply P, Q by s, t -> Lift x([s]*P), x([t]*Q) to full point -> Subtract [s]*P - [t]*Q -> Convert to PointX -> Create basis",
        |b| b.iter(|| {
            let [P_x, Q_x, ..] = params.two_torsion_basis.to_array();

            let P_x = params.starting_curve.xmul(&P_x, &s, s_bitsize);
            let Q_x = params.starting_curve.xmul(&Q_x, &s_inv, s_inv_bitsize);

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

                let P_x = params.starting_curve.xmul(&P_x, &s, s_bitsize);
                let Q_x = params.starting_curve.xmul(&Q_x, &s_inv, s_inv_bitsize);

                let (s_inv, _) = PokeFieldIBase::decode(&s_inv);
                let minus_s_inv = (PokeFieldIBase::MINUS_ONE * s_inv).encode();

                let PQ_x = params.starting_curve.ladder_biscalar(
                    &params.two_torsion_basis,
                    &s,
                    &minus_s_inv,
                    s_bitsize,
                    PokeFieldIBase::ENCODED_LENGTH,
                );

                BasisX::from_points(&P_x, &Q_x, &PQ_x)
            })
        },
    );
}

criterion_group!(benches, scalar_multiplication_then_basis_reconstruction);
criterion_main!(benches);
