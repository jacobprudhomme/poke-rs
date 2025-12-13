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

// POKE level III: p = 2^192 * 3^243 * 5^28 * 49 - 1
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

// POKE level V: p = 2^256 * 3^324 * 5^36 * 547 - 1
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

// INKE level I: p = 2^128 * 3^162 * 127 - 1
const INKE_I_MODULUS: [u64; 7] = [
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xbd94329674e2c7f6,
    0x0422d789525e827b,
    0xa7c33fd124a5e074,
    0xa8c4c5999e1b37e1,
    0x00000000000000d7,
];

fp2::define_fp2_from_modulus!(
    typename = InkeFieldI,
    base_typename = InkeFieldIBase,
    modulus = INKE_I_MODULUS,
);

// INKE level III: p = 2^192 * 3^243 * 5 * 7 - 1
const INKE_III_MODULUS: [u64; 10] = [
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0x6e6122aa036cb470,
    0x4a1cb212d489e150,
    0x1c5f1167f5c36a6c,
    0xa9f4e8afc92b3fe6,
    0x7975d4a56c357241,
    0x72e38a9838bafd79,
    0x000000000000004d,
];

fp2::define_fp2_from_modulus!(
    typename = InkeFieldIII,
    base_typename = InkeFieldIIIBase,
    modulus = INKE_III_MODULUS,
);

// INKE level V: p = 2^257 * 3^324 * 7^2 - 1
const INKE_V_MODULUS: [u64; 13] = [
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0xa83d19cb1c168101,
    0x5fb91904b4b3c20e,
    0x1a33717474b5ba03,
    0xb40334f01f15c342,
    0x73d4791f96fe9df3,
    0x84e53b06dea8956d,
    0x873beecd981293dd,
    0x96b6775d9e922f87,
    0x000000000000011a,
];

fp2::define_fp2_from_modulus!(
    typename = InkeFieldV,
    base_typename = InkeFieldVBase,
    modulus = INKE_V_MODULUS,
);

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
        fp2::define_fp2_tests!(PokeFieldIII, POKE_III_MODULUS, 6);
    }

    mod poke_v {
        use super::super::{POKE_V_MODULUS, PokeFieldV, PokeFieldVBase};

        fp2::define_fp_tests!(PokeFieldVBase);
        fp2::define_fp2_tests!(PokeFieldV, POKE_V_MODULUS, 5);
    }

    mod inke_i {
        use super::super::{INKE_I_MODULUS, InkeFieldI, InkeFieldIBase};

        fp2::define_fp_tests!(InkeFieldIBase);
        fp2::define_fp2_tests!(InkeFieldI, INKE_I_MODULUS, 2);
    }

    mod inke_iii {
        use super::super::{INKE_III_MODULUS, InkeFieldIII, InkeFieldIIIBase};

        fp2::define_fp_tests!(InkeFieldIIIBase);
        fp2::define_fp2_tests!(InkeFieldIII, INKE_III_MODULUS, 5);
    }

    mod inke_v {
        use super::super::{INKE_V_MODULUS, InkeFieldV, InkeFieldVBase};

        fp2::define_fp_tests!(InkeFieldVBase);
        fp2::define_fp2_tests!(InkeFieldV, INKE_V_MODULUS, 2);
    }
}
