#![allow(non_snake_case)]

use std::{io::Read as _, marker::PhantomData};

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::{
    elliptic::{basis::BasisX, curve::Curve},
    theta::elliptic_product::ProductPoint,
};
use ndarray::Array2;
use ndarray_rand::{RandomExt as _, rand::distributions::Uniform};
use num_bigint::{BigUint, RandBigInt as _};
use num_integer::gcd;
use sha3::{
    Shake256,
    digest::{ExtendableOutput as _, Update as _},
};

pub mod fields;
pub mod params;

pub const SUCCESS_RETVAL: u32 = u32::MAX;

pub struct PublicParams<Fp2: Fp2Trait> {
    pub starting_curve: Curve<Fp2>,
    pub two_torsion_exp: usize,
    pub three_torsion_exp: usize,
    pub five_torsion_exp: usize,
    pub two_torsion_basis: BasisX<Fp2>,
    pub three_torsion_basis: BasisX<Fp2>,
    pub five_torsion_basis: BasisX<Fp2>,
}

// FIXME: represent scalars as their LE byte arrays and bitsize. Removes external dependency on num-bigint
pub struct PrvKey<Fp2: Fp2Trait> {
    q: BigUint,
    alpha: BigUint,
    beta: BigUint,
    gamma: BigUint,
    delta: BigUint,
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

pub fn encrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    pub_key: &PubKey<Fp2>,
    message: &'a mut [u8],
) -> (Ciphertext<'a, Fp2>, u32) {
    // FIXME: where can I use vartime functions (i.e. operations on BigUint, gcd)? Where must things be constant-time?
    // TODO: add actual error handling

    let mut retval = SUCCESS_RETVAL;

    /* Sample scalars used for masking torsion points images or generating new kernels */

    let mut rng = rand::thread_rng();
    let ONE = BigUint::from(1u8);

    // The groups we will sample masking scalars from
    let two_torsion_order = BigUint::from(2u8).pow(
        pub_params
            .two_torsion_exp
            .try_into()
            .expect("Exponent of the 2-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let three_torsion_order = BigUint::from(3u8).pow(
        pub_params
            .three_torsion_exp
            .try_into()
            .expect("Exponent of the 3-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );
    let five_torsion_order = BigUint::from(5u8).pow(
        pub_params
            .five_torsion_exp
            .try_into()
            .expect("Exponent of the 5-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );

    // Sample scalar used to generate new kernels for sender's parallel isogenies
    let r = rng.gen_biguint_below(&three_torsion_order); // FIXME: what happens if this is 0?
    let r_bitsize =
        r.bits().try_into().expect("Size in bits of the scalar r is too big to fit in a usize (we do not ever expect this to happen)");
    let r = r.to_bytes_le();

    // Sample masking scalar for image of 2^a-torsion basis points on E_B and E_AB
    let mut omega = rng.gen_biguint_range(&ONE, &two_torsion_order);
    let mut omega_inv = omega.modinv(&two_torsion_order);
    while let None = omega_inv {
        println!("omega = {} is not invertible, retrying", omega);
        omega = rng.gen_biguint_range(&ONE, &two_torsion_order);
        omega_inv = omega.modinv(&two_torsion_order);
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
        Uniform::new(BigUint::ZERO, &five_torsion_order),
        &mut rng,
    );
    let mut det = (((&D[(0, 0)] * &D[(1, 1)]) % &five_torsion_order)
        + (&five_torsion_order - ((&D[(0, 1)] * &D[(1, 0)]) % &five_torsion_order)))
        % &five_torsion_order;
    let mut det_gcd = gcd(det.clone(), five_torsion_order.clone()); // TODO: look into a borrowing GCD function
    while det_gcd != ONE {
        println!("det(D) = {}, gcd(det(D), 5^c) = {}, retrying", det, det_gcd);
        D = Array2::random_using(
            (2, 2),
            Uniform::new(BigUint::ZERO, &five_torsion_order),
            &mut rng,
        );
        det = (((&D[(0, 0)] * &D[(1, 1)]) % &five_torsion_order)
            + (&five_torsion_order - ((&D[(0, 1)] * &D[(1, 0)]) % &five_torsion_order)))
            % &five_torsion_order;
        det_gcd = gcd(det.clone(), five_torsion_order.clone()); // TODO: look into a borrowing GCD function
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
    retval &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi to 2^a-torsion: {}",
        retval == SUCCESS_RETVAL,
    );

    let masked_xP_B = codomain_curve.xmul(P_B, &omega, omega_bitsize);
    let masked_xQ_B = codomain_curve.xmul(Q_B, &omega_inv, omega_inv_bitsize);

    let (masked_P_B, ok) = codomain_curve.lift_pointx(&masked_xP_B);
    retval &= ok;
    println!(
        "Successful execution after lifting x([omega]*P_B) to [omega]*P_B: {}",
        retval == SUCCESS_RETVAL,
    );
    let (masked_Q_B, ok) = codomain_curve.lift_pointx(&masked_xQ_B);
    retval &= ok;
    println!(
        "Successful execution after lifting x([1/omega]*Q_B) to [1/omega]*Q_B: {}",
        retval == SUCCESS_RETVAL,
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
    retval &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi to 5^c-torsion: {}",
        retval == SUCCESS_RETVAL,
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
    retval &= ok;
    println!(
        "Successful execution after lifting x(D*X_B) to D*X_B: {}",
        retval == SUCCESS_RETVAL,
    );
    let (masked_Y_B, ok) = codomain_curve_verif.lift_pointx(&masked_xY_B);
    retval &= ok;
    println!(
        "Successful execution after lifting x(D*Y_B) to D*Y_B: {}",
        retval == SUCCESS_RETVAL,
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
        SUCCESS_RETVAL,
    );

    // Apply sender's secret parallel isogeny to receiver's masked 2^a-torsion basis image points to obtain shared curve E_AB and pushforward basis image points (P_AB, Q_AB)
    let mut two_torsion_basis_EAB = pub_key.masked_two_torsion_basis_img.to_array();
    let (shared_end_curve, kernel_has_right_order) = pub_key.codomain_curve.three_isogeny_chain(
        &psi_prime_kernel,
        pub_params.three_torsion_exp,
        &mut two_torsion_basis_EAB,
    );
    let [xP_AB, xQ_AB, ..] = &two_torsion_basis_EAB;
    retval &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi' to 2^a-torsion: {}",
        retval == SUCCESS_RETVAL,
    );

    let masked_xP_AB = shared_end_curve.xmul(xP_AB, &omega, omega_bitsize);
    let masked_xQ_AB = shared_end_curve.xmul(xQ_AB, &omega_inv, omega_inv_bitsize);

    let (masked_P_AB, ok) = shared_end_curve.lift_pointx(&masked_xP_AB);
    retval &= ok;
    println!(
        "Successful execution after lifting x([omega]*P_AB) to [omega]*P_AB: {}",
        retval == SUCCESS_RETVAL,
    );
    let (masked_Q_AB, ok) = shared_end_curve.lift_pointx(&masked_xQ_AB);
    retval &= ok;
    println!(
        "Successful execution after lifting x([1/omega]*Q_AB) to [1/omega]*Q_AB: {}",
        retval == SUCCESS_RETVAL,
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
    retval &= kernel_has_right_order;
    println!(
        "Successful execution after applying psi' to 5^c-torsion: {}",
        retval == SUCCESS_RETVAL,
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
    retval &= ok;
    println!(
        "Successful execution after lifting x(D*X_AB) to D*X_AB: {}",
        retval == SUCCESS_RETVAL,
    );
    let (masked_Y_AB, ok) = shared_end_curve_verif.lift_pointx(&masked_xY_AB);
    retval &= ok;
    println!(
        "Successful execution after lifting x(D*Y_AB) to D*Y_AB: {}",
        retval == SUCCESS_RETVAL,
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
        SUCCESS_RETVAL,
    );

    let mut kdf = Shake256::default();
    let (X_AB, Y_AB) = shared_end_curve.lift_basis(&shared_secret);
    kdf.update(X_AB.to_string().as_bytes());
    kdf.update(Y_AB.to_string().as_bytes());
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
        retval,
    )
}

pub fn decrypt<'a, Fp2: Fp2Trait>(
    pub_params: &PublicParams<Fp2>,
    prv_key: &PrvKey<Fp2>,
    ciphertext: &Ciphertext<'a, Fp2>,
) -> (&'a [u8], u32) {
    let mut retval = SUCCESS_RETVAL;

    // The rings to work in to manipulate the scalars
    let two_torsion_order = BigUint::from(2u8).pow(
        pub_params
            .two_torsion_exp
            .try_into()
            .expect("Exponent of the 2-torsion subgroup is too big to fit in a u32 (we do not ever expect this to be the case)")
        );

    // Invert secret scalars, to neutralize their action on masked points we receive
    let alpha_inv = prv_key.alpha.modinv(&two_torsion_order);
    let Some(alpha_inv) = alpha_inv else {
        unreachable!();
    };
    let alpha_inv_bitsize =
        alpha_inv.bits().try_into().expect("Size in bits of the scalar 1/alpha is too big to fit in a usize (we do not ever expect this to happen)");
    let alpha_inv = alpha_inv.to_bytes_le();

    let beta_inv = prv_key.beta.modinv(&two_torsion_order);
    let Some(beta_inv) = beta_inv else {
        unreachable!();
    };
    let beta_inv_bitsize =
        beta_inv.bits().try_into().expect("Size in bits of the scalar 1/beta is too big to fit in a usize (we do not ever expect this to happen)");
    let beta_inv = beta_inv.to_bytes_le();

    // Neutralize action of our own secret scalars on the masked 2^a-torsion basis for E_AB
    let (P_AB, Q_AB) = ciphertext
        .shared_end_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EAB);
    let unmasked_P_AB = ciphertext
        .shared_end_curve
        .mul(&P_AB, &alpha_inv, alpha_inv_bitsize);
    let unmasked_Q_AB = ciphertext
        .shared_end_curve
        .mul(&Q_AB, &beta_inv, beta_inv_bitsize);

    // Compute kernel generators for our parallel 2D-isogeny Phi' (<([-q] P_B, P_AB'), ([-q] Q_B, Q_AB')>)
    let minus_q = &two_torsion_order - &prv_key.q;
    let minus_q_bitsize = minus_q.bits().try_into().expect("Size in bits of the scalar -q is too big to fit in a usize (we do not ever expect this to happen)");
    let minus_q = minus_q.to_bytes_le();

    let (P_B, Q_B) = ciphertext
        .codomain_curve
        .lift_basis(&ciphertext.masked_two_torsion_basis_EB);
    let generator_point1_B = ciphertext
        .codomain_curve
        .mul(&P_B, &minus_q, minus_q_bitsize);
    let generator_point2_B = ciphertext
        .codomain_curve
        .mul(&Q_B, &minus_q, minus_q_bitsize);

    let kernel_generator_point1 = ProductPoint::new(&generator_point1_B, &unmasked_P_AB);
    let kernel_generator_point2 = ProductPoint::new(&generator_point2_B, &unmasked_Q_AB);

    unimplemented!()
}
