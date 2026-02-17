#![allow(non_snake_case)]

use criterion::{Criterion, criterion_group, criterion_main};
use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, projective_point::Point};
use poke::{
    bn::BigNum,
    dlp::{
        solve_dlp_order_power_of_five, solve_dlp_order_power_of_five_explicit_subgroup,
        solve_dlp_small_prime_power_order,
    },
    params::{poke_i, poke_iii, poke_v},
    poke::PublicParams,
    rand::{sample_random_torsion_basis, sample_random_unit_mod_special_prime_power},
};

// Sample a basis of the 5^c-torsion, along with the Weil pairing on it,
// and a 3rd point we want to solve the discrete logarithm for
fn generate_fp2_power_of_five_torsion_basis_and_lc_point<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    pub_params: PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_23,
        NUM_WORDS_25,
        NUM_WORDS_35,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_2235,
        NUM_WORDS_2335,
        NUM_WORDS_2355,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
) -> (BasisX<Fp2>, Fp2, Point<Fp2>) {
    // The Weil pairing of a basis of the n-torsion subgroup will be the generator of an order-n subgroup of Fp^2
    let (U, V, eUV) = sample_random_torsion_basis(
        &pub_params.starting_curve,
        &[pub_params.five_torsion.base],
        &pub_params.five_torsion.order,
        &pub_params.five_torsion.cofactor,
    );
    let UV = pub_params.starting_curve.sub(&U, &V);

    let x = sample_random_unit_mod_special_prime_power(
        pub_params.five_torsion.base,
        &pub_params.five_torsion.order,
    );
    let y = sample_random_unit_mod_special_prime_power(
        pub_params.five_torsion.base,
        &pub_params.five_torsion.order,
    );

    let W1 = pub_params
        .starting_curve
        .mul(&U, &x.to_le_bytes(), x.nbits());
    let W2 = pub_params
        .starting_curve
        .mul(&V, &y.to_le_bytes(), y.nbits());
    let W = pub_params.starting_curve.add(&W1, &W2);

    (
        BasisX::from_points(&U.to_pointx(), &V.to_pointx(), &UV.to_pointx()),
        eUV,
        W,
    )
}

fn generic_method<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    pub_params: &PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_23,
        NUM_WORDS_25,
        NUM_WORDS_35,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_2235,
        NUM_WORDS_2335,
        NUM_WORDS_2355,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
    basis: &BasisX<Fp2>,
    basis_pairing: &Fp2,
    point: &Point<Fp2>,
) -> (BigNum<NUM_WORDS_5>, BigNum<NUM_WORDS_5>) {
    let (mut mU, V) = pub_params.starting_curve.lift_basis(basis);
    mU.set_neg();

    let WV = pub_params.starting_curve.sub(point, &V);
    let WmU = pub_params.starting_curve.sub(point, &mU);

    let covariant_pairing = pub_params.starting_curve.weil_pairing(
        &point.to_pointx().x(),
        &V.to_pointx().x(),
        &WV.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );
    let contravariant_pairing = pub_params.starting_curve.weil_pairing(
        &point.to_pointx().x(),
        &mU.to_pointx().x(),
        &WmU.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );

    let (x, _) = solve_dlp_small_prime_power_order(
        basis_pairing,
        &covariant_pairing,
        pub_params.five_torsion.base,
        pub_params.five_torsion.exp,
        &pub_params.five_torsion.p_adic_basis,
    );
    let (y, _) = solve_dlp_small_prime_power_order(
        basis_pairing,
        &contravariant_pairing,
        pub_params.five_torsion.base,
        pub_params.five_torsion.exp,
        &pub_params.five_torsion.p_adic_basis,
    );

    (x, y)
}

fn power_of_five_specialized_method<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    pub_params: &PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_23,
        NUM_WORDS_25,
        NUM_WORDS_35,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_2235,
        NUM_WORDS_2335,
        NUM_WORDS_2355,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
    basis: &BasisX<Fp2>,
    basis_pairing: &Fp2,
    point: &Point<Fp2>,
) -> (BigNum<NUM_WORDS_5>, BigNum<NUM_WORDS_5>) {
    let (mut mU, V) = pub_params.starting_curve.lift_basis(basis);
    mU.set_neg();

    let WV = pub_params.starting_curve.sub(point, &V);
    let WmU = pub_params.starting_curve.sub(point, &mU);

    let covariant_pairing = pub_params.starting_curve.weil_pairing(
        &point.to_pointx().x(),
        &V.to_pointx().x(),
        &WV.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );
    let contravariant_pairing = pub_params.starting_curve.weil_pairing(
        &point.to_pointx().x(),
        &mU.to_pointx().x(),
        &WmU.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );

    let (x, _) = solve_dlp_order_power_of_five(
        basis_pairing,
        &covariant_pairing,
        pub_params.five_torsion.exp,
        &pub_params.five_torsion.p_adic_basis,
    );
    let (y, _) = solve_dlp_order_power_of_five(
        basis_pairing,
        &contravariant_pairing,
        pub_params.five_torsion.exp,
        &pub_params.five_torsion.p_adic_basis,
    );

    (x, y)
}

fn power_of_five_specialized_method_explicit_subgroup<
    Fp2: Fp2Trait,
    const NUM_WORDS_2: usize,
    const NUM_WORDS_3: usize,
    const NUM_WORDS_5: usize,
    const NUM_WORDS_23: usize,
    const NUM_WORDS_25: usize,
    const NUM_WORDS_35: usize,
    const NUM_WORDS_P: usize,
    const NUM_WORDS_223: usize,
    const NUM_WORDS_2235: usize,
    const NUM_WORDS_2335: usize,
    const NUM_WORDS_2355: usize,
    const TWO_ADIC_BASIS_LEN: usize,
    const THREE_ADIC_BASIS_LEN: usize,
    const FIVE_ADIC_BASIS_LEN: usize,
>(
    pub_params: &PublicParams<
        Fp2,
        NUM_WORDS_2,
        NUM_WORDS_3,
        NUM_WORDS_5,
        NUM_WORDS_23,
        NUM_WORDS_25,
        NUM_WORDS_35,
        NUM_WORDS_P,
        NUM_WORDS_223,
        NUM_WORDS_2235,
        NUM_WORDS_2335,
        NUM_WORDS_2355,
        TWO_ADIC_BASIS_LEN,
        THREE_ADIC_BASIS_LEN,
        FIVE_ADIC_BASIS_LEN,
    >,
    basis: &BasisX<Fp2>,
    basis_pairing: &Fp2,
    point: &Point<Fp2>,
) -> (BigNum<NUM_WORDS_5>, BigNum<NUM_WORDS_5>) {
    let (U, V) = pub_params.starting_curve.lift_basis(basis);

    // Compute all points needed to explicitly construct the order-5 subgroup generated by e(U, V)^(5^(c-1))
    let minus_V = -V;
    let double_V = pub_params.starting_curve.double(&V);
    let minus_double_V = -double_V;
    let U2V = pub_params.starting_curve.sub(&U, &double_V);
    let Um2V = pub_params.starting_curve.sub(&U, &minus_double_V);
    let UmV = pub_params.starting_curve.sub(&U, &minus_V);

    // Compute all elements of the order-5 subgroup generated by e(U, V)^(5^(c-1))
    let generator_squared = pub_params.starting_curve.weil_pairing(
        &U.to_pointx().x(),
        &double_V.to_pointx().x(),
        &U2V.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );
    let generator_cubed = pub_params.starting_curve.weil_pairing(
        &U.to_pointx().x(),
        &minus_double_V.to_pointx().x(),
        &Um2V.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );
    let generator_to_four = pub_params.starting_curve.weil_pairing(
        &U.to_pointx().x(),
        &minus_V.to_pointx().x(),
        &UmV.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );
    let prime_order_subgroup = [
        Fp2::ONE,
        *basis_pairing,
        generator_squared,
        generator_cubed,
        generator_to_four,
    ];

    let mU = -U;
    let WV = pub_params.starting_curve.sub(point, &V);
    let WmU = pub_params.starting_curve.sub(point, &mU);

    let covariant_pairing = pub_params.starting_curve.weil_pairing(
        &point.to_pointx().x(),
        &V.to_pointx().x(),
        &WV.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );
    let contravariant_pairing = pub_params.starting_curve.weil_pairing(
        &point.to_pointx().x(),
        &mU.to_pointx().x(),
        &WmU.to_pointx().x(),
        &pub_params.five_torsion.order.to_le_bytes(),
        pub_params.five_torsion.order.nbits(),
    );

    let (x, _) = solve_dlp_order_power_of_five_explicit_subgroup(
        &prime_order_subgroup,
        &covariant_pairing,
        pub_params.five_torsion.exp,
        &pub_params.five_torsion.p_adic_basis,
    );
    let (y, _) = solve_dlp_order_power_of_five_explicit_subgroup(
        &prime_order_subgroup,
        &contravariant_pairing,
        pub_params.five_torsion.exp,
        &pub_params.five_torsion.p_adic_basis,
    );

    (x, y)
}

fn dlp(c: &mut Criterion) {
    let poke_i_params = poke_i::get_params();
    let poke_iii_params = poke_iii::get_params();
    let poke_v_params = poke_v::get_params();

    let (basis_i, basis_pairing_i, point_i) =
        generate_fp2_power_of_five_torsion_basis_and_lc_point(poke_i::get_params());
    let (basis_iii, basis_pairing_iii, point_iii) =
        generate_fp2_power_of_five_torsion_basis_and_lc_point(poke_iii::get_params());
    let (basis_v, basis_pairing_v, point_v) =
        generate_fp2_power_of_five_torsion_basis_and_lc_point(poke_v::get_params());

    let mut group_generic =
        c.benchmark_group("Discrete log in groups of prime order (= 5), using generic method");
    group_generic.bench_function("POKÉ level I", |b| {
        b.iter(|| generic_method(&poke_i_params, &basis_i, &basis_pairing_i, &point_i))
    });
    group_generic.bench_function("POKÉ level III", |b| {
        b.iter(|| generic_method(&poke_iii_params, &basis_iii, &basis_pairing_iii, &point_iii))
    });
    group_generic.bench_function("POKÉ level V", |b| {
        b.iter(|| generic_method(&poke_v_params, &basis_v, &basis_pairing_v, &point_v))
    });
    group_generic.finish();

    let mut group_specialized = c.benchmark_group(
        "Discrete log in groups of prime-power order (= 5^c), using specialized method",
    );
    group_specialized.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            power_of_five_specialized_method(&poke_i_params, &basis_i, &basis_pairing_i, &point_i)
        })
    });
    group_specialized.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            power_of_five_specialized_method(
                &poke_iii_params,
                &basis_iii,
                &basis_pairing_iii,
                &point_iii,
            )
        })
    });
    group_specialized.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            power_of_five_specialized_method(&poke_v_params, &basis_v, &basis_pairing_v, &point_v)
        })
    });
    group_specialized.finish();

    let mut group_specialized_explicit = c.benchmark_group(
        "Discrete log in groups of prime-power order (= 5^c), using specialized method where the order-5 subgroup is computed a priori",
    );
    group_specialized_explicit.bench_function("POKÉ level I", |b| {
        b.iter(|| {
            power_of_five_specialized_method_explicit_subgroup(
                &poke_i_params,
                &basis_i,
                &basis_pairing_i,
                &point_i,
            )
        })
    });
    group_specialized_explicit.bench_function("POKÉ level III", |b| {
        b.iter(|| {
            power_of_five_specialized_method_explicit_subgroup(
                &poke_iii_params,
                &basis_iii,
                &basis_pairing_iii,
                &point_iii,
            )
        })
    });
    group_specialized_explicit.bench_function("POKÉ level V", |b| {
        b.iter(|| {
            power_of_five_specialized_method_explicit_subgroup(
                &poke_v_params,
                &basis_v,
                &basis_pairing_v,
                &point_v,
            )
        })
    });
    group_specialized_explicit.finish();
}

criterion_group!(benches, dlp);
criterion_main!(benches);
