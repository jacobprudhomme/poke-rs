use core::marker::PhantomData;

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve, point::PointX, projective_point::Point},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
};
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _, XofReader as _},
};

use crate::{
    SUCCESS_RETVAL,
    bn::BigNum,
    dimtwo::{
        eval_2d_two_isogeny_chain_on_prime_power_torsion_basis,
        eval_2d_two_isogeny_chain_poke_separate_bases, generate_2d_isogeny_poke,
    },
    masking::{
        mask_basis_by_same_scalar, mask_basisx_by_diagonal_scalars,
        mask_basisx_by_diagonal_scalars_points_only, mask_basisx_by_same_scalar,
        mask_basisx_by_scalar_matrix, mask_basisx_by_scalar_matrix_pointx_only,
    },
    params::TorsionParams,
    rand::{
        sample_random_element_mod, sample_random_invertible_matrix_mod_prime_power,
        sample_random_secret_degree, sample_random_unit_mod_prime_power,
    },
};

pub struct PublicParams<
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
> {
    pub field_characteristic: BigNum<NUM_WORDS_P>,
    pub cofactor: BigNum<1>,
    pub starting_curve: Curve<Fp2>,
    pub effective_two_torsion_order: BigNum<NUM_WORDS_2>,
    pub effective_two_torsion_exp: usize,
    pub two_torsion: TorsionParams<NUM_WORDS_2, NUM_WORDS_35, TWO_ADIC_BASIS_LEN>,
    pub two_torsion_basis: BasisX<Fp2>,
    pub two_quaternion_basis: ([Point<Fp2>; 4], [Point<Fp2>; 4]),
    pub three_torsion: TorsionParams<NUM_WORDS_3, NUM_WORDS_25, THREE_ADIC_BASIS_LEN>,
    pub three_torsion_basis: BasisX<Fp2>,
    pub three_quaternion_basis: ([Point<Fp2>; 4], [Point<Fp2>; 4]),
    pub five_torsion: TorsionParams<NUM_WORDS_5, NUM_WORDS_23, FIVE_ADIC_BASIS_LEN>,
    pub five_torsion_basis: BasisX<Fp2>,
    pub full_torsion_order: BigNum<NUM_WORDS_P>,
    pub inv_three_order_mod_two_order: BigNum<NUM_WORDS_2>,
}

pub struct PrvKey<Fp2: Fp2Trait, const NUM_WORDS_2: usize, const NUM_WORDS_5: usize> {
    pub q: BigNum<NUM_WORDS_2>,
    pub alpha: BigNum<NUM_WORDS_2>,
    pub beta: BigNum<NUM_WORDS_2>,
    pub delta: BigNum<NUM_WORDS_5>,
    pub _field: PhantomData<Fp2>,
}

pub struct PubKey<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_img: BasisX<Fp2>,
    pub masked_three_torsion_basis_img: BasisX<Fp2>,
    pub masked_five_torsion_basis_img: BasisX<Fp2>,
}

pub struct Ciphertext<Fp2: Fp2Trait> {
    pub codomain_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EB: BasisX<Fp2>,
    pub masked_five_torsion_basis_EB: BasisX<Fp2>,
    pub shared_end_curve: Curve<Fp2>,
    pub masked_two_torsion_basis_EAB: BasisX<Fp2>,
    pub encrypted_message: Vec<u8>,
}

pub fn keygen<
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
) -> (PubKey<Fp2>, PrvKey<Fp2, NUM_WORDS_2, NUM_WORDS_5>, u32) {
    let mut retval = SUCCESS_RETVAL;

    // Keep generating elements until we find an invertible one
    // FIXME: I generate a q in [0, 2^(a-2) - 1] here. Should we be generating in [0, 2^(a-4)]
    // like in POKE-INKE? Or in [0, 2^(a-3) - 1] like in POKE-PKE?
    let (q, q_dual) = sample_random_secret_degree(
        &pub_params.effective_two_torsion_order,
        &[pub_params.three_torsion.base, pub_params.five_torsion.base],
    );

    let (domain, (P1P2, Q1Q2), ok) = generate_2d_isogeny_poke(pub_params, &q, &q_dual);
    retval &= ok;
    let codomain_curve = domain.curves().1;

    let ((P_A, Q_A), (R_A, S_A), (X_A, Y_A), ok) = eval_2d_two_isogeny_chain_poke_separate_bases(
        &domain,
        (&P1P2, &Q1Q2),
        pub_params.effective_two_torsion_exp,
        &q,
        &q_dual,
        (
            &pub_params.two_torsion,
            &pub_params.three_torsion,
            &pub_params.five_torsion,
        ),
        (
            &pub_params.two_torsion_basis,
            &pub_params.three_torsion_basis,
            &pub_params.five_torsion_basis,
        ),
    );
    retval &= ok;

    let two_torsion_basis_EA = BasisX::from_points(
        &P_A.to_pointx(),
        &Q_A.to_pointx(),
        &codomain_curve.sub(&P_A, &Q_A).to_pointx(),
    );
    let three_torsion_basis_EA = BasisX::from_points(
        &R_A.to_pointx(),
        &S_A.to_pointx(),
        &codomain_curve.sub(&R_A, &S_A).to_pointx(),
    );
    let five_torsion_basis_EA = BasisX::from_points(
        &X_A.to_pointx(),
        &Y_A.to_pointx(),
        &codomain_curve.sub(&X_A, &Y_A).to_pointx(),
    );

    let alpha = sample_random_unit_mod_prime_power(
        pub_params.two_torsion.base,
        &pub_params.two_torsion.order,
    );
    let beta = sample_random_unit_mod_prime_power(
        pub_params.two_torsion.base,
        &pub_params.two_torsion.order,
    );
    let gamma = sample_random_unit_mod_prime_power(
        pub_params.three_torsion.base,
        &pub_params.three_torsion.order,
    );
    let delta = sample_random_unit_mod_prime_power(
        pub_params.five_torsion.base,
        &pub_params.five_torsion.order,
    );

    let masked_two_torsion_basis_EA =
        mask_basisx_by_diagonal_scalars(&codomain_curve, &two_torsion_basis_EA, &alpha, &beta);
    let masked_three_torsion_basis_EA =
        mask_basisx_by_same_scalar(&codomain_curve, &three_torsion_basis_EA, &gamma);
    let masked_five_torsion_basis_EA =
        mask_basisx_by_same_scalar(&codomain_curve, &five_torsion_basis_EA, &delta);

    (
        PubKey {
            codomain_curve,
            masked_two_torsion_basis_img: masked_two_torsion_basis_EA,
            masked_three_torsion_basis_img: masked_three_torsion_basis_EA,
            masked_five_torsion_basis_img: masked_five_torsion_basis_EA,
        },
        PrvKey {
            q,
            alpha,
            beta,
            delta,
            _field: PhantomData,
        },
        retval,
    )
}

pub fn encrypt<
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
    pub_key: &PubKey<Fp2>,
    message: &[u8],
) -> (Ciphertext<Fp2>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    // FIXME: where can I use vartime functions (i.e. operations on BigUint, gcd)? Where must things be constant-time?

    let mut retval = SUCCESS_RETVAL;

    /* Sample scalars used for masking torsion points images or generating new kernels */

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = sample_random_element_mod(&pub_params.three_torsion.order);

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    // FIXME: Should this be full 2^a torsion, or effective 2^(a-2) torsion?
    let omega1 = sample_random_unit_mod_prime_power(
        pub_params.two_torsion.base,
        &pub_params.effective_two_torsion_order,
    );
    let omega2 = sample_random_unit_mod_prime_power(
        pub_params.two_torsion.base,
        &pub_params.effective_two_torsion_order,
    );

    // Sample masking matrix for image of 5^c-torsion basis points on E_B and E_AB
    let D = sample_random_invertible_matrix_mod_prime_power(5, &pub_params.five_torsion.order);

    /* Compute images of points, codomain curves through sender's secret parallel isogenies */

    // Compute kernel for sender's parallel isogenies psi (<R_0 + [r] S_0>) and psi' (<R_A + [r] S_A>)
    let psi_kernel = pub_params.starting_curve.three_point_ladder(
        &pub_params.three_torsion_basis,
        &r.to_le_bytes(),
        r.nbits(),
    );
    let psi_prime_kernel = pub_key.codomain_curve.three_point_ladder(
        &pub_key.masked_three_torsion_basis_img,
        &r.to_le_bytes(),
        r.nbits(),
    );

    // Apply sender's secret isogeny to 2^a- and 5^c-torsion bases to obtain their codomain curve E_B and basis image points (P_B, Q_B), (X_B, Y_B)
    let two_torsion_basis_EB = pub_params.two_torsion_basis.to_array();
    let five_torsion_basis_EB = pub_params.five_torsion_basis.to_array();
    let mut torsion_bases_EB = [PointX::INFINITY; 6];
    torsion_bases_EB[..3].copy_from_slice(&two_torsion_basis_EB);
    torsion_bases_EB[3..].copy_from_slice(&five_torsion_basis_EB);

    let (codomain_curve, kernel_has_right_order) = pub_params.starting_curve.three_isogeny_chain(
        &psi_kernel,
        pub_params.three_torsion.exp,
        &mut torsion_bases_EB,
    );
    retval &= kernel_has_right_order;

    let two_torsion_basis_EB = BasisX::from_slice(&torsion_bases_EB[..3]);
    let masked_two_torsion_basis_EB =
        mask_basisx_by_diagonal_scalars(&codomain_curve, &two_torsion_basis_EB, &omega1, &omega2);

    let five_torsion_basis_EB = BasisX::from_slice(&torsion_bases_EB[3..]);
    let masked_five_torsion_basis_EB =
        mask_basisx_by_scalar_matrix(&codomain_curve, &five_torsion_basis_EB, &D);

    // Apply sender's secret parallel isogeny to receiver's masked 2^a- and 5^c-torsion basis image points
    // to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB), (X_AB, Y_AB)
    let two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let five_torsion_basis_EAB = pub_key.masked_five_torsion_basis_img.to_array();
    let mut torsion_bases_EAB = [PointX::INFINITY; 6];
    torsion_bases_EAB[..3].copy_from_slice(&two_torsion_basis_EAB);
    torsion_bases_EAB[3..].copy_from_slice(&five_torsion_basis_EAB);

    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion.exp,
        &mut torsion_bases_EAB,
    );
    retval &= kernel_has_right_order;

    let two_torsion_basis_EAB = BasisX::from_slice(&torsion_bases_EAB[..3]);
    let masked_two_torsion_basis_EAB = mask_basisx_by_diagonal_scalars(
        &shared_end_curve,
        &two_torsion_basis_EAB,
        &omega1,
        &omega2,
    );

    let five_torsion_basis_EAB = BasisX::from_slice(&torsion_bases_EAB[3..]);
    let (masked_xX_AB, masked_xY_AB) =
        mask_basisx_by_scalar_matrix_pointx_only(&shared_end_curve, &five_torsion_basis_EAB, &D);

    let mut kdf = Shake256::default();
    kdf.update(&masked_xX_AB.x().encode());
    kdf.update(&masked_xY_AB.x().encode());
    let mut one_time_pad = kdf.finalize_xof();
    let mut encrypted_message = vec![0u8; message.len()];
    one_time_pad.read(&mut encrypted_message);
    for (encrypted_message_byte, message_byte) in encrypted_message.iter_mut().zip(message) {
        *encrypted_message_byte ^= message_byte;
    }

    let ct = Ciphertext {
        codomain_curve,
        masked_two_torsion_basis_EB,
        masked_five_torsion_basis_EB,
        shared_end_curve,
        masked_two_torsion_basis_EAB,
        encrypted_message,
    };

    (ct, retval)
}

pub fn decrypt<
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
    prv_key: &PrvKey<Fp2, NUM_WORDS_2, NUM_WORDS_5>,
    ciphertext: &Ciphertext<Fp2>,
) -> (Vec<u8>, u32)
where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut retval = SUCCESS_RETVAL;

    /* Construct kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>) */

    let alpha_q = prv_key
        .alpha
        .mul_mod_power_of_two(&prv_key.q, &pub_params.two_torsion.order);
    let beta_q = prv_key
        .beta
        .mul_mod_power_of_two(&prv_key.q, &pub_params.two_torsion.order);

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

    /* Derive shared secret from image points and use them to decrypt message */

    let q_dual = &pub_params.effective_two_torsion_order - &prv_key.q;
    let (_, five_torsion_basis_img_EAB, ok) =
        eval_2d_two_isogeny_chain_on_prime_power_torsion_basis(
            &domain,
            (&P1P2, &Q1Q2),
            pub_params.effective_two_torsion_exp,
            &prv_key.q,
            &q_dual,
            &pub_params.five_torsion,
            &ciphertext.masked_five_torsion_basis_EB,
        );
    retval &= ok;

    let (X_AB, Y_AB) = mask_basis_by_same_scalar(
        &ciphertext.shared_end_curve,
        &five_torsion_basis_img_EAB,
        &prv_key.delta,
    );

    // Undo one-time pad of message
    let mut kdf = Shake256::default();
    kdf.update(&X_AB.to_pointx().x().encode());
    kdf.update(&Y_AB.to_pointx().x().encode());
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
