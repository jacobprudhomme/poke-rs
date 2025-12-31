use std::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _, XofReader as _},
};

use crate::{
    SUCCESS_RETVAL,
    bn::BigNum,
    masking::{mask_basisx_by_diagonal_scalars, mask_basisx_by_diagonal_scalars_points_only},
    rand::{sample_random_element_mod, sample_random_unit_mod_prime_power},
};

pub struct PublicParams<Fp2: Fp2Trait, const NUM_WORDS_2: usize, const NUM_WORDS_3: usize> {
    pub starting_curve: Curve<Fp2>,
    pub full_two_torsion_order: BigNum<NUM_WORDS_2>,
    pub full_two_torsion_exp: usize,
    pub effective_two_torsion_order: BigNum<NUM_WORDS_2>,
    pub effective_two_torsion_exp: usize,
    pub three_torsion_order: BigNum<NUM_WORDS_3>,
    pub three_torsion_exp: usize,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
}

pub struct PrvKey<Fp2: Fp2Trait, const NUM_WORDS_2: usize> {
    pub q: BigNum<NUM_WORDS_2>,
    pub alpha: BigNum<NUM_WORDS_2>,
    pub beta: BigNum<NUM_WORDS_2>,
    pub _field: PhantomData<Fp2>,
}

pub struct PubKey<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_img: BasisX<Fp2>,
    pub masked_three_torsion_basis_img: BasisX<Fp2>,
    pub intermediate_curve: Curve<Fp2>,
    pub masked_three_torsion_basis_img_intermediate: BasisX<Fp2>,
}

pub struct Ciphertext<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EB: BasisX<Fp2>,
    pub shared_end_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EAB: BasisX<Fp2>,
    pub encrypted_message: Vec<u8>,
}

pub fn encrypt<Fp2: Fp2Trait, const NUM_WORDS_2: usize, const NUM_WORDS_3: usize>(
    pub_params: &PublicParams<Fp2, NUM_WORDS_2, NUM_WORDS_3>,
    pub_key: &PubKey<Fp2>,
    message: &[u8],
) -> (Ciphertext<Fp2>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    /* Sample scalars used for masking torsion points images or generating new kernels */

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = sample_random_element_mod(&pub_params.three_torsion_order);

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    let omega1 = sample_random_unit_mod_prime_power(2, &pub_params.effective_two_torsion_order);
    let omega2 = sample_random_unit_mod_prime_power(2, &pub_params.effective_two_torsion_order);

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>), psi' (<R_A + [r] S_A>) and psi'' (<R_A1 + [r] S_A1>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r.to_le_bytes(),
        r.nbits(),
    );
    let psi_prime_kernel = pub_key.intermediate_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img_intermediate,
        &r.to_le_bytes(),
        r.nbits(),
    );
    let psi_dblprime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.to_le_bytes(),
        r.nbits(),
    );

    // Apply sender's secret isogeny to 2^a-torsion basis to obtain their codomain curve E_B and basis image points (P_B, Q_B)
    let mut two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EB,
    );
    retval &= kernel_has_right_order;

    let two_torsion_basis_EB = BasisX::from_slice(&two_torsion_basis_EB);
    let masked_two_torsion_basis_EB =
        mask_basisx_by_diagonal_scalars(&codomain_curve, &two_torsion_basis_EB, &omega1, &omega2);

    // Apply sender's secret parallel isogeny to receiver's masked 2^a-torsion basis image points to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB)
    let mut two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_dblprime_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EAB,
    );
    retval &= kernel_has_right_order;

    let two_torsion_basis_EAB = BasisX::from_slice(&two_torsion_basis_EAB);
    let masked_two_torsion_basis_EAB = mask_basisx_by_diagonal_scalars(
        &shared_end_curve,
        &two_torsion_basis_EAB,
        &omega1,
        &omega2,
    );

    // Compute codomain of sender's secret intermediate parallel isogeny to obtain shared secret curve
    let (secret_curve, _) = pub_key.intermediate_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut [],
    );

    // Compute shared secret from j-invariant of shared secret curve and encrypt message
    let mut kdf = Shake256::default();
    kdf.update(&secret_curve.j_invariant().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut encrypted_message = vec![0u8; message.len()];
    one_time_pad.read(&mut encrypted_message);
    for (encrypted_message_byte, message_byte) in encrypted_message.iter_mut().zip(message) {
        *encrypted_message_byte ^= message_byte;
    }

    let ct = Ciphertext {
        codomain_curve,
        masked_two_torsion_basis_EB,
        shared_end_curve,
        masked_two_torsion_basis_EAB,
        encrypted_message,
    };

    (ct, retval)
}

pub fn decrypt<Fp2: Fp2Trait, const NUM_WORDS_2: usize, const NUM_WORDS_3: usize>(
    pub_params: &PublicParams<Fp2, NUM_WORDS_2, NUM_WORDS_3>,
    prv_key: &PrvKey<Fp2, NUM_WORDS_2>,
    ciphertext: &Ciphertext<Fp2>,
) -> (Vec<u8>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    // Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let alpha_q = prv_key
        .alpha
        .mul_mod_power_of_two(&prv_key.q, &pub_params.full_two_torsion_order);
    let beta_q = prv_key
        .beta
        .mul_mod_power_of_two(&prv_key.q, &pub_params.full_two_torsion_order);

    let (scaled_P_B, scaled_Q_B) = mask_basisx_by_diagonal_scalars_points_only(
        &ciphertext.codomain_curve,
        &ciphertext.masked_two_torsion_basis_EB,
        &alpha_q,
        &beta_q,
    );

    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);

    let domain = EllipticProduct::new(&ciphertext.codomain_curve, &ciphertext.shared_end_curve);
    let P1P2 = ProductPoint::new(&scaled_P_B, &P_AB);
    let Q1Q2 = ProductPoint::new(&scaled_Q_B, &Q_AB);

    // Compute codomain curve pair of Phi', which contains the shared secret curve
    let (aux_curves, _, ok) =
        domain.elliptic_product_isogeny(&P1P2, &Q1Q2, pub_params.effective_two_torsion_exp, &[]);
    retval &= ok;
    let secret_curve = aux_curves.curves().0;

    // Undo one-time pad of message
    let mut kdf = Shake256::default();
    kdf.update(&secret_curve.j_invariant().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut message = vec![0u8; ciphertext.encrypted_message.len()];
    one_time_pad.read(&mut message);
    for (message_byte, encrypted_message_byte) in
        message.iter_mut().zip(&ciphertext.encrypted_message)
    {
        *message_byte ^= encrypted_message_byte;
    }

    (message, retval)
}
