#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use poke::{
    FAILURE_RETVAL,
    inke::{decrypt, encrypt, keygen},
    params,
};
use rand::RngCore;

fn inke(c: &mut Criterion) {
    let params_i = params::inke_i::get_params();
    let params_iii = params::inke_iii::get_params();
    let params_v = params::inke_v::get_params();

    let (mut pub_key_i, mut prv_key_i, mut ok) = keygen(&params_i);
    while ok == FAILURE_RETVAL {
        (pub_key_i, prv_key_i, ok) = keygen(&params_i);
    }
    let (mut pub_key_iii, mut prv_key_iii, mut ok) = keygen(&params_iii);
    while ok == FAILURE_RETVAL {
        (pub_key_iii, prv_key_iii, ok) = keygen(&params_iii);
    }
    let (mut pub_key_v, mut prv_key_v, mut ok) = keygen(&params_v);
    while ok == FAILURE_RETVAL {
        (pub_key_v, prv_key_v, ok) = keygen(&params_v);
    }

    let mut rng = rand::rng();
    let mut message = [0; 128];
    rng.fill_bytes(&mut message);

    let (ct_i, _) = encrypt(&params_i, &pub_key_i, &message);
    let (ct_iii, _) = encrypt(&params_iii, &pub_key_iii, &message);
    let (ct_v, _) = encrypt(&params_v, &pub_key_v, &message);

    let mut keygen_group = c.benchmark_group("Key Generation");
    keygen_group.bench_function("INKE level I", |b| b.iter(|| keygen(&params_i)));
    keygen_group.bench_function("INKE level III", |b| b.iter(|| keygen(&params_iii)));
    keygen_group.bench_function("INKE level V", |b| b.iter(|| keygen(&params_v)));
    keygen_group.finish();

    let mut encryption_group = c.benchmark_group("Encryption");
    encryption_group.bench_function("INKE level I", |b| {
        b.iter(|| encrypt(&params_i, &pub_key_i, &message))
    });
    encryption_group.bench_function("INKE level III", |b| {
        b.iter(|| encrypt(&params_iii, &pub_key_iii, &message))
    });
    encryption_group.bench_function("INKE level V", |b| {
        b.iter(|| encrypt(&params_v, &pub_key_v, &message))
    });
    encryption_group.finish();

    let mut decryption_group = c.benchmark_group("Decryption");
    decryption_group.bench_function("INKE level I", |b| {
        b.iter(|| decrypt(&params_i, &prv_key_i, &ct_i))
    });
    decryption_group.bench_function("INKE level III", |b| {
        b.iter(|| decrypt(&params_iii, &prv_key_iii, &ct_iii))
    });
    decryption_group.bench_function("INKE level V", |b| {
        b.iter(|| decrypt(&params_v, &prv_key_v, &ct_v))
    });
    decryption_group.finish();
}

criterion_group!(benches, inke);
criterion_main!(benches);
