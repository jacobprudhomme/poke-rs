#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use poke::{
    example_keypairs, params,
    poke::{decrypt, encrypt},
};
use rand::RngCore;

fn poke(c: &mut Criterion) {
    let params = params::poke_i::get_params();
    let pub_key = example_keypairs::poke_i::get_pub_key();
    let prv_key = example_keypairs::poke_i::get_prv_key();

    let mut rng = rand::rng();
    let mut message = [0; 128];
    rng.fill_bytes(&mut message);

    let mut encryption_group = c.benchmark_group("Encryption");
    encryption_group.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            encrypt(&params, &pub_key, &message);
        })
    });
    encryption_group.finish();

    let (ct, _) = encrypt(&params, &pub_key, &message);

    let mut decryption_group = c.benchmark_group("Decryption");
    decryption_group.bench_function("POKÉ level I", |b| {
        b.iter(|| decrypt(&params, &prv_key, &ct));
    });
    decryption_group.finish();
}

criterion_group!(benches, poke);
criterion_main!(benches);
