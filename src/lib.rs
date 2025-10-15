#![allow(non_snake_case)]

use std::{io::Read as _, marker::PhantomData};

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, curve::Curve};
use ndarray::Array2;
use ndarray_rand::{RandomExt as _, rand::distributions::Uniform};
use num_bigint::{BigUint, RandBigInt as _};
use num_integer::gcd;
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _},
};

// POKE level I: p = 2^129 * 3^164 * 5^18 - 1
const POKE_I_MODULUS: [u64; 7] = [
    0xffffffffffffffff,
    0xffffffffffffffff,
    0x3d346b3e65f69451,
    0x7ef3ecff193099d0,
    0x56ff93faead91477,
    0xc6124673c50d17a5,
    0x00006a0bf4180690,
];

fp2::define_fp2_from_modulus!(
    typename = PokeFieldI,
    base_typename = PokeFieldIBase,
    modulus = POKE_I_MODULUS,
);

// POKE level II: p = 2^192 * 3^243 * 5^28 * 49 - 1
const POKE_III_MODULUS: [u64; 11] = [
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xab4c1ec4e9a4421a,
    0xa1a751e0ff03064a,
    0x5c5381a82432b77b,
    0x74f54cc513a36773,
    0x152ef0c01f75ccd4,
    0xa53054622a07450c,
    0xf81dcb46fd3f8b4d,
    0x00000000000000da,
];

fp2::define_fp2_from_modulus!(
    typename = PokeFieldIII,
    base_typename = PokeFieldIIIBase,
    modulus = POKE_III_MODULUS,
);

// POKE level III: p = 2^256 * 3^324 * 5^36 * 547 - 1
const POKE_V_MODULUS: [u64; 14] = [
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xe8334b6ad7209ce2,
    0x0912211ae1688d39,
    0x926e7814cc1dd2be,
    0x370d3afd3477d13d,
    0x2a5efc5fc51c347f,
    0x44282a1040b00581,
    0x61d48d4313219044,
    0x462f78afb014b0f7,
    0x507c1448c8213857,
    0x000000004a2a22b7,
];

fp2::define_fp2_from_modulus!(
    typename = PokeFieldV,
    base_typename = PokeFieldVBase,
    modulus = POKE_V_MODULUS,
);

pub struct PublicParams<Fp2: Fp2Trait> {
    starting_curve: Curve<Fp2>,
    two_torsion_exp: usize,
    three_torsion_exp: usize,
    five_torsion_exp: usize,
    two_torsion_basis: BasisX<Fp2>,
    three_torsion_basis: BasisX<Fp2>,
    five_torsion_basis: BasisX<Fp2>,
}

pub struct PrvKey<'a, Fp2: Fp2Trait> {
    q: usize,
    alpha: &'a [u8],
    beta: &'a [u8],
    gamma: &'a [u8],
    delta: &'a [u8],
    _field: PhantomData<Fp2>,
}

pub struct PubKey<Fp2: Fp2Trait> {
    codomain_curve: Curve<Fp2>,
    masked_two_torsion_basis_img: BasisX<Fp2>,
    masked_three_torsion_basis_img: BasisX<Fp2>,
    masked_five_torsion_basis_img: BasisX<Fp2>,
}

pub struct Ciphertext<'a, Fp2: Fp2Trait> {
    codomain_curve: Curve<Fp2>,
    masked_two_torsion_basis_img: BasisX<Fp2>,
    masked_five_torsion_basis_img: BasisX<Fp2>,
    shared_codomain_curve: Curve<Fp2>,
    masked_two_torsion_basis_pushfwd_img: BasisX<Fp2>,
    encrypted_message: &'a [u8],
}

pub fn encrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    pub_key: &PubKey<Fp2>,
    message: &'a mut [u8],
) -> Ciphertext<'a, Fp2> {
    /* Sample scalars used for masking torsion points images or generating new kernels */

    let mut rng = rand::thread_rng();
    let ONE = BigUint::from(1u8);

    // The subgroups we will sample from
    let Z_two_torsion_order = BigUint::from(2u8).pow(
        pub_params
            .two_torsion_exp
            .try_into()
            .expect("Exponent of the 2-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let Z_three_torsion_order = BigUint::from(3u8).pow(
        pub_params
            .three_torsion_exp
            .try_into()
            .expect("Exponent of the 3-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let Z_five_torsion_order = BigUint::from(5u8).pow(
        pub_params
            .five_torsion_exp
            .try_into()
            .expect("Exponent of the 5-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = rng.gen_biguint_below(&Z_three_torsion_order); // FIXME: what happens if this is 0?
    let r_bitsize =
        r.bits().try_into().expect("Size in bits of the scalar r is too big to fit in a usize (we do not ever expect this to happen)");
    let r = r.to_bytes_le();

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    let mut omega = rng.gen_biguint_range(&ONE, &Z_two_torsion_order);
    let mut omega_inv = omega.modinv(&Z_two_torsion_order);
    while let None = omega_inv {
        println!("omega = {} is not invertible, retrying", omega);
        omega = rng.gen_biguint_range(&ONE, &Z_two_torsion_order);
        omega_inv = omega.modinv(&Z_two_torsion_order);
    }
    println!();
    let Some(omega_inv) = omega_inv else {
        unreachable!();
    };
    let omega_bitsize =
        omega.bits().try_into().expect("Size in bits of the scalar omega is too big to fit in a usize (we do not ever expect this to happen)");
    let omega_inv_bitsize =
        omega_inv.bits().try_into().expect("Size in bits of the scalar 1/omega is too big to fit in a usize (we do not ever expect this to happen)");
    let omega = omega.to_bytes_le();
    let omega_inv = omega_inv.to_bytes_le();

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    // FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
    let mut D = Array2::random_using(
        (2, 2),
        Uniform::new(BigUint::ZERO, &Z_five_torsion_order),
        &mut rng,
    );
    let mut det = (((&D[(0, 0)] * &D[(1, 1)]) % &Z_five_torsion_order)
        + (&Z_five_torsion_order - ((&D[(0, 1)] * &D[(1, 0)]) % &Z_five_torsion_order)))
        % &Z_five_torsion_order;
    let mut det_gcd = gcd(det.clone(), Z_five_torsion_order.clone()); // TODO: look into a borrowing GCD function
    while det_gcd != ONE {
        println!("det(D) = {}, gcd(det(D), 5^c) = {}, retrying", det, det_gcd);
        D = Array2::random_using(
            (2, 2),
            Uniform::new(BigUint::ZERO, &Z_five_torsion_order),
            &mut rng,
        );
        det = (((&D[(0, 0)] * &D[(1, 1)]) % &Z_five_torsion_order)
            + (&Z_five_torsion_order - ((&D[(0, 1)] * &D[(1, 0)]) % &Z_five_torsion_order)))
            % &Z_five_torsion_order;
        det_gcd = gcd(det.clone(), Z_five_torsion_order.clone()); // TODO: look into a borrowing GCD function
    }
    println!();
    let D_bitsize = D.map(|elem| {
        TryInto::<usize>::try_into(elem.bits())
            .expect("Size in bits of the scalar is too big to fit in a usize (we do not ever expect this to happen)")
    });
    let D = D.map(|elem| elem.to_bytes_le());

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>) and psi' (<R_A + [r] S_A>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r,
        r_bitsize,
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r,
        r_bitsize,
    );

    // Apply sender's secret isogeny to 2^a-torsion basis to obtain their codomain curve E_B and basis image points (P_B, Q_B)
    // FIXME: there must be a more straightforward way to operate on individual points in an x-only basis than to
    // lift it to a standard basis -> obtain new P - Q point -> create new x-only basis from x-coordinates
    let mut two_torsion_basis_img = pub_params.two_torsion_basis.to_array();
    let (codomain_curve, _) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_img,
    );
    let [P_img, Q_img, ..] = &two_torsion_basis_img;
    let masked_P_img_x = codomain_curve.xmul(P_img, &omega, omega_bitsize);
    let masked_Q_img_x = codomain_curve.xmul(Q_img, &omega_inv, omega_inv_bitsize);
    let (masked_P_img, _) = codomain_curve.lift_pointx(&masked_P_img_x);
    let (masked_Q_img, _) = codomain_curve.lift_pointx(&masked_Q_img_x);
    let masked_PQ_img = codomain_curve.sub(&masked_P_img, &masked_Q_img);
    let masked_two_torsion_basis_img =
        BasisX::from_points(&masked_P_img_x, &masked_Q_img_x, &masked_PQ_img.to_pointx());

    // Apply sender's secret isogeny to 5^c-torsion basis to obtain basis image points (X_B, Y_B)
    // FIXME: there must be a more straightforward way to operate on individual points in an x-only basis than to
    // lift it to a standard basis -> obtain new P - Q point -> create new x-only basis from x-coordinates
    let mut five_torsion_basis_img = pub_params.five_torsion_basis.to_array();
    let (codomain_curve_verif, _) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut five_torsion_basis_img,
    );
    let five_torsion_basis_img = BasisX::from_slice(&five_torsion_basis_img);
    let masked_X_img_x = codomain_curve_verif.ladder_biscalar(
        &five_torsion_basis_img,
        &D[(0, 0)],
        &D[(0, 1)],
        D_bitsize[(0, 0)],
        D_bitsize[(0, 1)],
    );
    let (masked_X_img, _) = codomain_curve_verif.lift_pointx(&masked_X_img_x);
    let masked_Y_img_x = codomain_curve_verif.ladder_biscalar(
        &five_torsion_basis_img,
        &D[(1, 0)],
        &D[(1, 1)],
        D_bitsize[(1, 0)],
        D_bitsize[(1, 1)],
    );
    let (masked_Y_img, _) = codomain_curve_verif.lift_pointx(&masked_Y_img_x);
    let masked_XY_img = codomain_curve_verif.sub(&masked_X_img, &masked_Y_img);
    let masked_five_torsion_basis_img =
        BasisX::from_points(&masked_X_img_x, &masked_Y_img_x, &masked_XY_img.to_pointx());

    println!("j-invariant for sender's codomain curve:");
    println!("{}", codomain_curve.j_invariant());
    println!("{}\n", codomain_curve_verif.j_invariant());
    assert_eq!(
        codomain_curve
            .j_invariant()
            .equals(&codomain_curve_verif.j_invariant()),
        0xffffffff
    );

    // Apply sender's secret parallel isogeny to receiver's masked 2^a-torsion basis image points to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB)
    // FIXME: there must be a more straightforward way to operate on individual points in an x-only basis than to
    // lift it to a standard basis -> obtain new P - Q point -> create new x-only basis from x-coordinates
    let mut two_torsion_basis_pushfwd_img =
        pub_key.masked_two_torsion_basis_img.to_array();
    let (shared_codomain_curve, _) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_pushfwd_img,
    );
    let [P_AB, Q_AB, ..] = &two_torsion_basis_pushfwd_img;
    let masked_P_AB_x = shared_codomain_curve.xmul(P_AB, &omega, omega_bitsize);
    let masked_Q_AB_x = shared_codomain_curve.xmul(Q_AB, &omega_inv, omega_inv_bitsize);
    let (masked_P_AB, _) = shared_codomain_curve.lift_pointx(&masked_P_AB_x);
    let (masked_Q_AB, _) = shared_codomain_curve.lift_pointx(&masked_Q_AB_x);
    let masked_PQ_AB = shared_codomain_curve.sub(&masked_P_AB, &masked_Q_AB);
    let masked_two_torsion_basis_pushfwd_img =
        BasisX::from_points(&masked_P_AB_x, &masked_Q_AB_x, &masked_PQ_AB.to_pointx());

    // Apply sender's secret parallel isogeny to receiver's masked 5^c-torsion basis image points to obtain shared secret (X_AB, Y_AB)
    // FIXME: there must be a more straightforward way to operate on individual points in an x-only basis than to
    // lift it to a standard basis -> obtain new P - Q point -> create new x-only basis from x-coordinates
    let mut five_torsion_basis_pushfwd_img =
        pub_key.masked_five_torsion_basis_img.to_array();
    let (shared_codomain_curve_verif, _) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut five_torsion_basis_pushfwd_img,
    );
    let shared_secret =
        BasisX::from_slice(&five_torsion_basis_pushfwd_img);
    let masked_X_AB_x = shared_codomain_curve_verif.ladder_biscalar(
        &shared_secret,
        &D[(0, 0)],
        &D[(0, 1)],
        D_bitsize[(0, 0)],
        D_bitsize[(0, 1)],
    );
    let (masked_X_AB, _) = shared_codomain_curve_verif.lift_pointx(&masked_X_AB_x);
    let masked_Y_AB_x = shared_codomain_curve_verif.ladder_biscalar(
        &shared_secret,
        &D[(1, 0)],
        &D[(1, 1)],
        D_bitsize[(1, 0)],
        D_bitsize[(1, 1)],
    );
    let (masked_Y_AB, _) = shared_codomain_curve_verif.lift_pointx(&masked_Y_AB_x);
    let masked_XY_AB = shared_codomain_curve_verif.sub(&masked_X_AB, &masked_Y_AB);
    let shared_secret =
        BasisX::from_points(&masked_X_AB_x, &masked_Y_AB_x, &masked_XY_AB.to_pointx());

    println!("j-invariant for the shared end curve:");
    println!("{}", shared_codomain_curve.j_invariant());
    println!("{}\n", shared_codomain_curve_verif.j_invariant());
    assert_eq!(
        shared_codomain_curve
            .j_invariant()
            .equals(&shared_codomain_curve_verif.j_invariant()),
        0xffffffff
    );

    let mut kdf = Shake256::default();
    kdf.update(shared_secret.P.to_string().as_bytes());
    kdf.update(shared_secret.Q.to_string().as_bytes());
    let mut one_time_pad = kdf.finalize_xof();
    let mut buffer = vec![0u8; message.len()];
    let Ok(_) = one_time_pad.read(&mut buffer) else {
        panic!("Could not read bytes from KDF");
    };
    for (message_byte, one_time_pad_byte) in message.iter_mut().zip(buffer) {
        *message_byte ^= one_time_pad_byte;
    }

    Ciphertext {
        codomain_curve,
        masked_two_torsion_basis_img,
        masked_five_torsion_basis_img,
        shared_codomain_curve,
        masked_two_torsion_basis_pushfwd_img,
        encrypted_message: message,
    }
}

#[cfg(test)]
mod tests {
    mod poke_i {
        use super::super::{POKE_I_MODULUS, PokeFieldI, PokeFieldIBase};

        fp2::define_fp_tests!(PokeFieldIBase);
        fp2::define_fp2_tests!(PokeFieldI, POKE_I_MODULUS, 5);
    }

    mod poke_iii {
        use super::super::{POKE_III_MODULUS, PokeFieldIII, PokeFieldIIIBase};

        fp2::define_fp_tests!(PokeFieldIIIBase);
        fp2::define_fp2_tests!(PokeFieldIII, POKE_III_MODULUS, 5);
    }

    mod poke_v {
        use super::super::{POKE_V_MODULUS, PokeFieldV, PokeFieldVBase};

        fp2::define_fp_tests!(PokeFieldVBase);
        fp2::define_fp2_tests!(PokeFieldV, POKE_V_MODULUS, 5);
    }
}
