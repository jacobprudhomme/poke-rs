pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

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

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert!(true);
    }

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
