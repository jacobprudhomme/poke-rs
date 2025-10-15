#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]
#![recursion_limit = "256"]

use std::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, curve::Curve};
use ndarray::arr2;
use rand::RngCore as _;

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
    pushthrough_codomain_curve: Curve<Fp2>,
    masked_two_torsion_basis_pushthrough_img: BasisX<Fp2>,
    encrypted_message: &'a [u8],
}

pub fn encrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    pub_key: &PubKey<Fp2>,
    message: &'a mut [u8],
) -> Ciphertext<'a, Fp2>
    where [(); Fp2::ENCODED_LENGTH]: Sized
{
    /* Sample scalars used for masking torsion points images or generating new kernels */

    let mut rng = rand::thread_rng();

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    // FIXME: implement proper sampling of this value
    let mut r = [0u8; Fp2::ENCODED_LENGTH]; // FIXME: find exact bit size for elements in Z_(3^b)
    rng.fill_bytes(&mut r);

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    // FIXME: implement proper sampling of this value (have no way of inverting it without passing back to a bigint)
    let mut omega = [0u8; Fp2::ENCODED_LENGTH]; // FIXME: find exact bit size for elements in Z_(2^a), and make sure only invertible elements are allowed
    rng.fill_bytes(&mut omega);

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    // FIXME: implement proper sampling of this value (find algorithms to generate uniformly random determinant-1 matrices in SL_2(Z_(5^c)))
    let mut D = arr2(&[[[0u8; Fp2::ENCODED_LENGTH]; 2]; 2]); // FIXME: find exact bit size for elements in Z_(5^c), and make sure determinant is 1
    rng.fill_bytes(&mut D[(0, 0)]);
    rng.fill_bytes(&mut D[(0, 1)]);
    rng.fill_bytes(&mut D[(1, 0)]);
    rng.fill_bytes(&mut D[(1, 1)]);

    unimplemented!()
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
