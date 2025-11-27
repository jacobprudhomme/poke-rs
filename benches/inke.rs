#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use poke::{example_keypairs, inke_decrypt, inke_encrypt, params};
use rand::RngCore;

fn inke(c: &mut Criterion) {
    let params = params::inke_i::get_params();
    let pub_key = example_keypairs::inke_i::get_pub_key();
    let prv_key = example_keypairs::inke_i::get_prv_key();

    let mut rng = rand::rng();
    let mut message = [0; 128];
    rng.fill_bytes(&mut message);

    let mut encryption_group = c.benchmark_group("Encryption");
    encryption_group.bench_function("INKE level I", |b| {
        b.iter(|| {
            inke_encrypt(&params, &pub_key, &message);
        })
    });
    encryption_group.finish();

    let (ct, _) = inke_encrypt(&params, &pub_key, &message);

    let mut decryption_group = c.benchmark_group("Decryption");
    decryption_group.bench_function("INKE level I", |b| {
        b.iter(|| inke_decrypt(&params, &prv_key, &ct));
    });
    decryption_group.finish();
}

criterion_group!(benches, inke);
criterion_main!(benches);
