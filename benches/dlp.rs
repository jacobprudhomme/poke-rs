#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use fp2::traits::{Fp as _, Fp2 as Fp2Trait};
use poke::{
    bn::BigNum,
    dlp::{
        solve_dlp_order_five, solve_dlp_order_power_of_five, solve_dlp_small_prime_order,
        solve_dlp_small_prime_power_order,
    },
    params::{poke_i, poke_iii, poke_v},
    poke::PublicParams,
    rand::sample_random_torsion_basis,
};

// Sample a generator of a subgroup of Fp^2 of order 5^exp
fn generate_fp2_power_of_five_subgroup_generator<Fp2: Fp2Trait>(
    pub_params: PublicParams<Fp2>,
    exp: usize,
) -> Fp2 {
    let torsion_order = BigNum::from_prime_power(5, exp);
    let cofactor_order = &pub_params.five_torsion_cofactor
        * &BigNum::from_prime_power(5, pub_params.five_torsion_exp - exp);

    // The Weil pairing of a basis of the n-torsion subgroup will be the generator of an order-n subgroup of Fp^2
    let (_, _, eUV) =
        sample_random_torsion_basis(&pub_params.starting_curve, &torsion_order, &cofactor_order);

    eUV
}

fn dlp_prime_order(c: &mut Criterion) {
    const PRIME: usize = 5;
    const EXP: usize = 1;
    const LOG: [u8; 1] = [3];
    const LOG_BITSIZE: usize = 2;

    let gen_i = generate_fp2_power_of_five_subgroup_generator(poke_i::get_params(), EXP);
    let val_i = gen_i.pow(&LOG, LOG_BITSIZE);
    let gen_iii = generate_fp2_power_of_five_subgroup_generator(poke_iii::get_params(), EXP);
    let val_iii = gen_iii.pow(&LOG, LOG_BITSIZE);
    let gen_v = generate_fp2_power_of_five_subgroup_generator(poke_v::get_params(), EXP);
    let val_v = gen_v.pow(&LOG, LOG_BITSIZE);

    let mut group_generic =
        c.benchmark_group("Discrete log in groups of prime order (= 5), using generic method");
    group_generic.bench_function("POKÉ level I", |b| {
        b.iter(|| solve_dlp_small_prime_order(&gen_i, &val_i, PRIME))
    });
    group_generic.bench_function("POKÉ level III", |b| {
        b.iter(|| solve_dlp_small_prime_order(&gen_iii, &val_iii, PRIME))
    });
    group_generic.bench_function("POKÉ level V", |b| {
        b.iter(|| solve_dlp_small_prime_order(&gen_v, &val_v, PRIME))
    });
    group_generic.finish();

    let mut group_specialized =
        c.benchmark_group("Discrete log in groups of prime order (= 5), using specialized method");
    group_specialized.bench_function("POKÉ level I", |b| {
        b.iter(|| solve_dlp_order_five(&gen_i, &val_i))
    });
    group_specialized.bench_function("POKÉ level III", |b| {
        b.iter(|| solve_dlp_order_five(&gen_iii, &val_iii))
    });
    group_specialized.bench_function("POKÉ level V", |b| {
        b.iter(|| solve_dlp_order_five(&gen_v, &val_v))
    });
    group_specialized.finish();
}

fn dlp_prime_power_order(c: &mut Criterion) {
    const PRIME: usize = 5;
    const EXP: usize = 13;
    const LOG: [u8; 4] = [21, 205, 91, 7]; // 123456789
    const LOG_BITSIZE: usize = 27;

    let gen_i = generate_fp2_power_of_five_subgroup_generator(poke_i::get_params(), EXP);
    let val_i = gen_i.pow(&LOG, LOG_BITSIZE);
    let gen_iii = generate_fp2_power_of_five_subgroup_generator(poke_iii::get_params(), EXP);
    let val_iii = gen_iii.pow(&LOG, LOG_BITSIZE);
    let gen_v = generate_fp2_power_of_five_subgroup_generator(poke_v::get_params(), EXP);
    let val_v = gen_v.pow(&LOG, LOG_BITSIZE);

    let mut group_generic = c.benchmark_group(
        "Discrete log in groups of prime-power order (= 5^c), using generic method",
    );
    group_generic.bench_function("POKÉ level I", |b| {
        b.iter(|| solve_dlp_small_prime_power_order(&gen_i, &val_i, PRIME, EXP))
    });
    group_generic.bench_function("POKÉ level III", |b| {
        b.iter(|| solve_dlp_small_prime_power_order(&gen_iii, &val_iii, PRIME, EXP))
    });
    group_generic.bench_function("POKÉ level V", |b| {
        b.iter(|| solve_dlp_small_prime_power_order(&gen_v, &val_v, PRIME, EXP))
    });
    group_generic.finish();

    let mut group_specialized = c.benchmark_group(
        "Discrete log in groups of prime-power order (= 5^c), using specialized method",
    );
    group_specialized.bench_function("POKÉ level I", |b| {
        b.iter(|| solve_dlp_order_power_of_five(&gen_i, &val_i, EXP))
    });
    group_specialized.bench_function("POKÉ level III", |b| {
        b.iter(|| solve_dlp_order_power_of_five(&gen_iii, &val_iii, EXP))
    });
    group_specialized.bench_function("POKÉ level V", |b| {
        b.iter(|| solve_dlp_order_power_of_five(&gen_v, &val_v, EXP))
    });
    group_specialized.finish();
}

criterion_group!(benches, dlp_prime_order, dlp_prime_power_order);
criterion_main!(benches);
