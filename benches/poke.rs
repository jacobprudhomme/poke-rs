#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use poke::{
    example_keypairs, params,
    poke::{decrypt, encrypt},
};
use rand::RngCore;

fn poke(c: &mut Criterion) {
    let params_i = params::poke_i::get_params();
    let pub_key_i = example_keypairs::poke_i::get_pub_key();
    let prv_key_i = example_keypairs::poke_i::get_prv_key();
    let params_v = params::poke_v::get_params();
    let pub_key_v = example_keypairs::poke_v::get_pub_key();
    let prv_key_v = example_keypairs::poke_v::get_prv_key();

    let mut rng = rand::rng();
    let mut message = [0; 128];
    rng.fill_bytes(&mut message);

    let mut encryption_group = c.benchmark_group("Encryption");
    encryption_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            encrypt(&params_i, &pub_key_i, &message);
        })
    });
    encryption_group.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            encrypt(&params_v, &pub_key_v, &message);
        })
    });
    encryption_group.finish();

    let (ct_i, _) = encrypt(&params_i, &pub_key_i, &message);
    let (ct_v, _) = encrypt(&params_v, &pub_key_v, &message);

    let mut decryption_group = c.benchmark_group("Decryption");
    decryption_group.bench_function("POKÉ level I", |b| {
        b.iter(|| decrypt(&params_i, &prv_key_i, &ct_i));
    });
    decryption_group.bench_function("POKÉ level V", |b| {
        b.iter(|| decrypt(&params_v, &prv_key_v, &ct_v));
    });
    decryption_group.finish();
}

criterion_group!(benches, poke);
criterion_main!(benches);
