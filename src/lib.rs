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
    pub starting_curve: Curve<Fp2>,
    pub two_torsion_exp: usize,
    pub three_torsion_exp: usize,
    pub five_torsion_exp: usize,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
    pub five_torsion_basis: BasisX<Fp2>,
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
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_img: BasisX<Fp2>,
    pub masked_three_torsion_basis_img: BasisX<Fp2>,
    pub masked_five_torsion_basis_img: BasisX<Fp2>,
}

pub struct Ciphertext<'a, Fp2: Fp2Trait> {
    codomain_curve: Curve<Fp2>,
    masked_two_torsion_basis_EB: BasisX<Fp2>,
    masked_five_torsion_basis_EB: BasisX<Fp2>,
    shared_end_curve: Curve<Fp2>,
    masked_two_torsion_basis_EAB: BasisX<Fp2>,
    pub encrypted_message: &'a [u8],
}

pub mod poke_i {
    use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX};

    use crate::{PokeFieldI, PublicParams};

    // Construct basis points for the 2^a-torsion on E_0
    const P_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        225, 204, 48, 233, 181, 190, 18, 222, 247, 38, 59, 93, 252, 209, 65, 62, 195, 253, 222, 58,
        179, 18, 119, 130, 98, 196, 148, 139, 59, 204, 93, 73, 22, 7, 63, 63, 184, 164, 108, 255,
        205, 79, 133, 20, 182, 27, 46, 205, 220, 82, 131, 215, 39, 28,
    ];
    const P_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        156, 1, 55, 113, 211, 191, 79, 224, 97, 54, 107, 37, 254, 167, 210, 138, 199, 125, 108,
        159, 62, 27, 61, 12, 176, 93, 127, 206, 236, 40, 77, 235, 18, 81, 163, 191, 61, 216, 30,
        105, 141, 244, 112, 38, 122, 199, 207, 251, 158, 170, 70, 187, 238, 73,
    ];
    const Q_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        225, 204, 48, 233, 181, 190, 18, 222, 247, 38, 59, 93, 252, 209, 65, 62, 195, 253, 222, 58,
        179, 18, 119, 130, 98, 196, 148, 139, 59, 204, 93, 73, 22, 7, 63, 63, 184, 164, 108, 255,
        205, 79, 133, 20, 182, 27, 46, 205, 220, 82, 131, 215, 39, 28,
    ];
    const Q_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        99, 254, 200, 142, 44, 64, 176, 31, 158, 201, 148, 218, 1, 88, 45, 117, 138, 22, 138, 198,
        255, 79, 247, 48, 32, 60, 177, 74, 18, 196, 166, 147, 100, 195, 53, 43, 189, 187, 224, 237,
        23, 35, 156, 158, 249, 126, 66, 202, 241, 91, 209, 56, 29, 32,
    ];
    const PQ_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        84, 5, 39, 47, 28, 26, 72, 53, 222, 25, 244, 169, 18, 138, 123, 250, 227, 91, 135, 191,
        182, 168, 208, 156, 231, 66, 10, 171, 57, 90, 207, 9, 222, 195, 240, 102, 7, 222, 148, 122,
        208, 175, 249, 130, 55, 245, 12, 92, 175, 174, 252, 231, 208, 15,
    ];
    const PQ_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [0; PokeFieldI::ENCODED_LENGTH / 2];
    const P_X: PokeFieldI = PokeFieldI::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: PokeFieldI = PokeFieldI::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Construct basis points for the 3^b-torsion on E_0
    const R_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        33, 100, 192, 219, 131, 122, 237, 66, 95, 99, 60, 177, 230, 250, 51, 190, 104, 113, 44,
        242, 139, 87, 147, 181, 249, 53, 197, 220, 252, 127, 88, 234, 23, 241, 221, 97, 160, 52,
        102, 44, 37, 165, 139, 203, 245, 120, 204, 216, 248, 102, 186, 121, 47, 14,
    ];
    const R_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        242, 78, 54, 243, 244, 158, 40, 209, 213, 36, 144, 132, 126, 115, 146, 252, 115, 95, 79,
        49, 121, 90, 228, 120, 114, 82, 233, 129, 214, 22, 113, 22, 116, 81, 115, 222, 238, 180,
        157, 29, 159, 205, 134, 216, 253, 65, 214, 79, 148, 149, 147, 24, 195, 7,
    ];
    const S_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        33, 100, 192, 219, 131, 122, 237, 66, 95, 99, 60, 177, 230, 250, 51, 190, 104, 113, 44,
        242, 139, 87, 147, 181, 249, 53, 197, 220, 252, 127, 88, 234, 23, 241, 221, 97, 160, 52,
        102, 44, 37, 165, 139, 203, 245, 120, 204, 216, 248, 102, 186, 121, 47, 14,
    ];
    const S_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        13, 177, 201, 12, 11, 97, 215, 46, 42, 219, 111, 123, 129, 140, 109, 3, 222, 52, 167, 52,
        197, 16, 80, 196, 93, 71, 71, 151, 40, 214, 130, 104, 3, 195, 101, 12, 12, 223, 97, 57, 6,
        74, 134, 236, 117, 4, 60, 118, 252, 112, 132, 219, 72, 98,
    ];
    const RS_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        214, 198, 30, 107, 72, 179, 73, 69, 35, 50, 116, 38, 27, 143, 85, 161, 55, 176, 109, 176,
        64, 247, 227, 127, 52, 115, 253, 72, 217, 177, 78, 213, 224, 192, 75, 192, 253, 45, 130,
        177, 170, 220, 184, 89, 185, 137, 120, 89, 231, 163, 80, 255, 92, 95,
    ];
    const RS_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [0; PokeFieldI::ENCODED_LENGTH / 2];
    const R_X: PokeFieldI = PokeFieldI::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: PokeFieldI = PokeFieldI::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: PokeFieldI = PokeFieldI::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    // Construct basis points for the 5^c-torsion on E_0
    const X_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        142, 141, 154, 162, 86, 251, 208, 110, 83, 81, 167, 239, 99, 27, 248, 99, 176, 209, 50, 79,
        95, 226, 187, 103, 115, 94, 168, 239, 128, 125, 222, 127, 12, 58, 148, 85, 96, 16, 38, 236,
        30, 216, 153, 163, 196, 201, 222, 27, 117, 237, 189, 56, 217, 95,
    ];
    const X_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        211, 24, 39, 150, 242, 9, 42, 168, 202, 224, 82, 61, 102, 182, 231, 124, 213, 107, 144, 72,
        84, 15, 181, 210, 65, 156, 234, 60, 141, 56, 253, 222, 254, 41, 3, 136, 237, 101, 182, 89,
        189, 117, 17, 158, 8, 209, 192, 197, 185, 255, 80, 133, 107, 67,
    ];
    const Y_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        142, 141, 154, 162, 86, 251, 208, 110, 83, 81, 167, 239, 99, 27, 248, 99, 176, 209, 50, 79,
        95, 226, 187, 103, 115, 94, 168, 239, 128, 125, 222, 127, 12, 58, 148, 85, 96, 16, 38, 236,
        30, 216, 153, 163, 196, 201, 222, 27, 117, 237, 189, 56, 217, 95,
    ];
    const Y_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        44, 231, 216, 105, 13, 246, 213, 87, 53, 31, 173, 194, 153, 73, 24, 131, 124, 40, 102, 29,
        234, 91, 127, 106, 142, 253, 69, 220, 113, 180, 246, 159, 120, 234, 213, 98, 13, 46, 73,
        253, 231, 161, 251, 38, 107, 117, 81, 0, 215, 6, 199, 110, 160, 38,
    ];
    const XY_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        7, 127, 166, 170, 67, 68, 236, 216, 13, 61, 70, 222, 190, 115, 147, 244, 207, 140, 116,
        141, 195, 61, 63, 202, 239, 236, 93, 15, 92, 242, 111, 151, 53, 67, 144, 196, 218, 77, 91,
        160, 138, 225, 199, 32, 138, 40, 48, 229, 231, 107, 234, 178, 109, 61,
    ];
    const XY_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [0; PokeFieldI::ENCODED_LENGTH / 2];
    const X_X: PokeFieldI = PokeFieldI::const_decode_no_check(&X_X_RE, &X_X_IM);
    const Y_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Y_X_RE, &Y_X_IM);
    const XY_X: PokeFieldI = PokeFieldI::const_decode_no_check(&XY_X_RE, &XY_X_IM);

    pub fn create_poke_i_params() -> PublicParams<PokeFieldI> {
        let two_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&P_X),
            &PointX::from_x_coord(&Q_X),
            &PointX::from_x_coord(&PQ_X),
        );
        let three_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&R_X),
            &PointX::from_x_coord(&S_X),
            &PointX::from_x_coord(&RS_X),
        );
        let five_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&X_X),
            &PointX::from_x_coord(&Y_X),
            &PointX::from_x_coord(&XY_X),
        );

        PublicParams {
            starting_curve: Curve::new(&PokeFieldI::from_i32(6)),
            two_torsion_exp: 129,
            three_torsion_exp: 164,
            five_torsion_exp: 18,
            two_torsion_basis,
            three_torsion_basis,
            five_torsion_basis,
        }
    }
}

pub fn encrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    pub_key: &PubKey<Fp2>,
    message: &'a mut [u8],
) -> (Ciphertext<'a, Fp2>, u32) {
    // FIXME: where can I use vartime functions (i.e. operations on BigUint, gcd)? Where must things be constant-time?
    // TODO: add actual error handling

    let mut retcode = 0xffffffff;

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
    let mut two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EB,
    );
    let [P_B, Q_B, ..] = &two_torsion_basis_EB;
    retcode &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi to 2^a-torsion: {}",
        retcode == 0xffffffff,
    );

    let masked_xP_B = codomain_curve.xmul(P_B, &omega, omega_bitsize);
    let masked_xQ_B = codomain_curve.xmul(Q_B, &omega_inv, omega_inv_bitsize);

    let (masked_P_B, ok) = codomain_curve.lift_pointx(&masked_xP_B);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x([omega]*P_B) to [omega]*P_B: {}",
        retcode == 0xffffffff,
    );
    let (masked_Q_B, ok) = codomain_curve.lift_pointx(&masked_xQ_B);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x([1/omega]*Q_B) to [1/omega]*Q_B: {}",
        retcode == 0xffffffff,
    );

    let masked_PQ_B = codomain_curve.sub(&masked_P_B, &masked_Q_B);

    let masked_two_torsion_basis_EB =
        BasisX::from_points(&masked_xP_B, &masked_xQ_B, &masked_PQ_B.to_pointx());

    // Apply sender's secret isogeny to 5^c-torsion basis to obtain basis image points (X_B, Y_B)
    let mut five_torsion_basis_EB = pub_params.five_torsion_basis.to_array();
    let (codomain_curve_verif, kernel_has_right_order) =
        pub_params.starting_curve.three_isogeny_chain(
            &psi_kernel,
            pub_params.three_torsion_exp,
            &mut five_torsion_basis_EB,
        );
    let five_torsion_basis_EB = BasisX::from_slice(&five_torsion_basis_EB);
    retcode &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi to 5^c-torsion: {}",
        retcode == 0xffffffff,
    );

    let masked_xX_B = codomain_curve_verif.ladder_biscalar(
        &five_torsion_basis_EB,
        &D[(0, 0)],
        &D[(0, 1)],
        D_bitsize[(0, 0)],
        D_bitsize[(0, 1)],
    );
    let masked_xY_B = codomain_curve_verif.ladder_biscalar(
        &five_torsion_basis_EB,
        &D[(1, 0)],
        &D[(1, 1)],
        D_bitsize[(1, 0)],
        D_bitsize[(1, 1)],
    );

    let (masked_X_B, ok) = codomain_curve_verif.lift_pointx(&masked_xX_B);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x(D*X_B) to D*X_B: {}",
        retcode == 0xffffffff,
    );
    let (masked_Y_B, ok) = codomain_curve_verif.lift_pointx(&masked_xY_B);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x(D*Y_B) to D*Y_B: {}",
        retcode == 0xffffffff,
    );

    let masked_XY_B = codomain_curve_verif.sub(&masked_X_B, &masked_Y_B);

    let masked_five_torsion_basis_EB =
        BasisX::from_points(&masked_xX_B, &masked_xY_B, &masked_XY_B.to_pointx());

    println!("j-invariant for sender's codomain curve:");
    println!("{}", codomain_curve.j_invariant());
    println!("{}\n", codomain_curve_verif.j_invariant());
    assert_eq!(
        codomain_curve
            .j_invariant()
            .equals(&codomain_curve_verif.j_invariant()),
        0xffffffff,
    );

    // Apply sender's secret parallel isogeny to receiver's masked 2^a-torsion basis image points to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB)
    let mut two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EAB,
    );
    let [xP_AB, xQ_AB, ..] = &two_torsion_basis_EAB;
    retcode &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi' to 2^a-torsion: {}",
        retcode == 0xffffffff,
    );

    let masked_xP_AB = shared_end_curve.xmul(xP_AB, &omega, omega_bitsize);
    let masked_xQ_AB = shared_end_curve.xmul(xQ_AB, &omega_inv, omega_inv_bitsize);

    let (masked_P_AB, ok) = shared_end_curve.lift_pointx(&masked_xP_AB);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x([omega]*P_AB) to [omega]*P_AB: {}",
        retcode == 0xffffffff,
    );
    let (masked_Q_AB, ok) = shared_end_curve.lift_pointx(&masked_xQ_AB);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x([1/omega]*Q_AB) to [1/omega]*Q_AB: {}",
        retcode == 0xffffffff,
    );

    let masked_PQ_AB = shared_end_curve.sub(&masked_P_AB, &masked_Q_AB);

    let masked_two_torsion_basis_EAB =
        BasisX::from_points(&masked_xP_AB, &masked_xQ_AB, &masked_PQ_AB.to_pointx());

    // Apply sender's secret parallel isogeny to receiver's masked 5^c-torsion basis image points to obtain shared secret (X_AB, Y_AB)
    let mut five_torsion_basis_EAB = pub_key.masked_five_torsion_basis_img.to_array();
    let (shared_end_curve_verif, kernel_has_right_order) =
        pub_key.codomain_curve.three_isogeny_chain(
            &psi_prime_kernel,
            pub_params.three_torsion_exp,
            &mut five_torsion_basis_EAB,
        );
    let shared_secret = BasisX::from_slice(&five_torsion_basis_EAB);
    retcode &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi' to 5^c-torsion: {}",
        retcode == 0xffffffff,
    );

    let masked_xX_AB = shared_end_curve_verif.ladder_biscalar(
        &shared_secret,
        &D[(0, 0)],
        &D[(0, 1)],
        D_bitsize[(0, 0)],
        D_bitsize[(0, 1)],
    );
    let masked_xY_AB = shared_end_curve_verif.ladder_biscalar(
        &shared_secret,
        &D[(1, 0)],
        &D[(1, 1)],
        D_bitsize[(1, 0)],
        D_bitsize[(1, 1)],
    );

    let (masked_X_AB, ok) = shared_end_curve_verif.lift_pointx(&masked_xX_AB);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x(D*X_AB) to D*X_AB: {}",
        retcode == 0xffffffff,
    );
    let (masked_Y_AB, ok) = shared_end_curve_verif.lift_pointx(&masked_xY_AB);
    retcode &= ok;
    println!(
        "Successful execution after applying lifting x(D*Y_AB) to D*Y_AB: {}",
        retcode == 0xffffffff,
    );

    let masked_XY_AB = shared_end_curve_verif.sub(&masked_X_AB, &masked_Y_AB);

    let shared_secret =
        BasisX::from_points(&masked_xX_AB, &masked_xY_AB, &masked_XY_AB.to_pointx());

    println!("j-invariant for the shared end curve:");
    println!("{}", shared_end_curve.j_invariant());
    println!("{}\n", shared_end_curve_verif.j_invariant());
    assert_eq!(
        shared_end_curve
            .j_invariant()
            .equals(&shared_end_curve_verif.j_invariant()),
        0xffffffff,
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

    (
        Ciphertext {
            codomain_curve,
            masked_two_torsion_basis_EB,
            masked_five_torsion_basis_EB,
            shared_end_curve,
            masked_two_torsion_basis_EAB,
            encrypted_message: message,
        },
        retcode,
    )
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
