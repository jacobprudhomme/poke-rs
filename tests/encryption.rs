#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]

use fp2::traits::Fp2 as Fp2Trait;
use poke::{
    FAILURE_RETVAL, SUCCESS_RETVAL,
    params::{poke_i, poke_iii, poke_v},
    poke::{PubKey, PublicParams, encrypt, keygen},
};
use rstest::rstest;
use rstest_reuse::{apply, template};

const MESSAGE: &'static [u8; 13] = b"Hello, world!";

fn pub_key<
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
    params: PublicParams<
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
) -> PubKey<Fp2> {
    let (mut pub_key, _, mut ok) = keygen(&params);
    while ok == FAILURE_RETVAL {
        (pub_key, _, ok) = keygen(&params);
    }

    pub_key
}

#[template]
#[rstest]
#[case::poke_level_i(poke_i::get_params(), pub_key(poke_i::get_params()))]
#[case::poke_level_iii(poke_iii::get_params(), pub_key(poke_iii::get_params()))]
#[case::poke_level_v(poke_v::get_params(), pub_key(poke_v::get_params()))]
fn encryption_test_data<Fp2: Fp2Trait>(
    #[case] params: PublicParams<Fp2>,
    #[case] pub_key: PubKey<Fp2>,
) {
}

#[apply(encryption_test_data)]
fn encryption_passes<
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
    params: PublicParams<
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
    pub_key: PubKey<Fp2>,
) where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut message = *MESSAGE;
    let message = message.as_mut_slice();
    let (_, ok) = encrypt(&params, &pub_key, message);

    assert_eq!(ok, SUCCESS_RETVAL);
}

#[apply(encryption_test_data)]
fn message_and_ciphertext_have_same_length<
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
    params: PublicParams<
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
    pub_key: PubKey<Fp2>,
) where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut message = *MESSAGE;
    let message = message.as_mut_slice();
    let (ct, _) = encrypt(&params, &pub_key, message);

    assert_eq!(ct.encrypted_message.len(), MESSAGE.len());
}

#[apply(encryption_test_data)]
fn masked_image_points_have_correct_order<
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
    params: PublicParams<
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
    pub_key: PubKey<Fp2>,
) where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut message = *MESSAGE;
    let message = message.as_mut_slice();
    let (ct, _) = encrypt(&params, &pub_key, message);

    // Ensure [2^a] * (P_B, Q_B, P_B - Q_B) = O, and [2^a - 1] * (P_B, Q_B, P_B - Q_B) != O
    let two_torsion_basis_times_one_less_than_order = ct
        .codomain_curve
        .basis_double_iter(&ct.masked_two_torsion_basis_EB, params.two_torsion.exp - 1);
    assert_eq!(
        two_torsion_basis_times_one_less_than_order.P.is_zero(),
        FAILURE_RETVAL
    );
    assert_eq!(
        two_torsion_basis_times_one_less_than_order.Q.is_zero(),
        FAILURE_RETVAL
    );
    assert_eq!(
        two_torsion_basis_times_one_less_than_order.PQ.is_zero(),
        FAILURE_RETVAL
    );

    let two_torsion_basis_times_order = ct
        .codomain_curve
        .basis_double_iter(&ct.masked_two_torsion_basis_EB, params.two_torsion.exp);
    assert_eq!(two_torsion_basis_times_order.P.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.Q.is_zero(), SUCCESS_RETVAL);
    assert_eq!(two_torsion_basis_times_order.PQ.is_zero(), SUCCESS_RETVAL);

    // Ensure [2^a] * (P_AB, Q_AB, P_AB - Q_AB) = O, and [2^a - 1] * (P_AB, Q_AB, P_AB - Q_AB) != O
    let xP_AB_times_one_less_than_order = ct.shared_end_curve.xmul(
        &ct.masked_two_torsion_basis_EAB.P,
        &params.two_torsion.reduced_order.to_le_bytes(),
        params.two_torsion.reduced_order.nbits(),
    );
    let xQ_AB_times_one_less_than_order = ct.shared_end_curve.xmul(
        &ct.masked_two_torsion_basis_EAB.Q,
        &params.two_torsion.reduced_order.to_le_bytes(),
        params.two_torsion.reduced_order.nbits(),
    );
    let xPQ_AB_times_one_less_than_order = ct.shared_end_curve.xmul(
        &ct.masked_two_torsion_basis_EAB.PQ,
        &params.two_torsion.reduced_order.to_le_bytes(),
        params.two_torsion.reduced_order.nbits(),
    );
    assert_eq!(xP_AB_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xQ_AB_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xPQ_AB_times_one_less_than_order.is_zero(), FAILURE_RETVAL);

    let xP_AB_times_order = ct.shared_end_curve.xmul(
        &ct.masked_two_torsion_basis_EAB.P,
        &params.two_torsion.order.to_le_bytes(),
        params.two_torsion.order.nbits(),
    );
    let xQ_AB_times_order = ct.shared_end_curve.xmul(
        &ct.masked_two_torsion_basis_EAB.Q,
        &params.two_torsion.order.to_le_bytes(),
        params.two_torsion.order.nbits(),
    );
    let xPQ_AB_times_order = ct.shared_end_curve.xmul(
        &ct.masked_two_torsion_basis_EAB.PQ,
        &params.two_torsion.order.to_le_bytes(),
        params.two_torsion.order.nbits(),
    );
    assert_eq!(xP_AB_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xQ_AB_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xPQ_AB_times_order.is_zero(), SUCCESS_RETVAL);

    // Ensure [5^c] * (X_B, Y_B, X_B - Y_B) = O and [5^c - 1] * (X_B, Y_B, X_B - Y_B) != O
    let xX_times_one_less_than_order = ct.codomain_curve.xmul(
        &ct.masked_five_torsion_basis_EB.P,
        &params.five_torsion.reduced_order.to_le_bytes(),
        params.five_torsion.reduced_order.nbits(),
    );
    let xY_times_one_less_than_order = ct.codomain_curve.xmul(
        &ct.masked_five_torsion_basis_EB.Q,
        &params.five_torsion.reduced_order.to_le_bytes(),
        params.five_torsion.reduced_order.nbits(),
    );
    let xXY_times_one_less_than_order = ct.codomain_curve.xmul(
        &ct.masked_five_torsion_basis_EB.PQ,
        &params.five_torsion.reduced_order.to_le_bytes(),
        params.five_torsion.reduced_order.nbits(),
    );
    assert_eq!(xX_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xY_times_one_less_than_order.is_zero(), FAILURE_RETVAL);
    assert_eq!(xXY_times_one_less_than_order.is_zero(), FAILURE_RETVAL);

    let xX_times_order = ct.codomain_curve.xmul(
        &ct.masked_five_torsion_basis_EB.P,
        &params.five_torsion.order.to_le_bytes(),
        params.five_torsion.order.nbits(),
    );
    let xY_times_order = ct.codomain_curve.xmul(
        &ct.masked_five_torsion_basis_EB.Q,
        &params.five_torsion.order.to_le_bytes(),
        params.five_torsion.order.nbits(),
    );
    let xXY_times_order = ct.codomain_curve.xmul(
        &ct.masked_five_torsion_basis_EB.PQ,
        &params.five_torsion.order.to_le_bytes(),
        params.five_torsion.order.nbits(),
    );
    assert_eq!(xX_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xY_times_order.is_zero(), SUCCESS_RETVAL);
    assert_eq!(xXY_times_order.is_zero(), SUCCESS_RETVAL);
}

#[apply(encryption_test_data)]
fn masked_image_points_are_on_curve<
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
    params: PublicParams<
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
    pub_key: PubKey<Fp2>,
) where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut message = *MESSAGE;
    let message = message.as_mut_slice();
    let (ct, _) = encrypt(&params, &pub_key, message);

    // Check (P_B, Q_B, P_B - Q_B), a basis of the 2^a-torsion on E_B
    assert_eq!(
        ct.codomain_curve
            .is_on_curve(&ct.masked_two_torsion_basis_EB.P.x()),
        SUCCESS_RETVAL,
        "P_B is not on the curve E_B",
    );
    assert_eq!(
        ct.codomain_curve
            .is_on_curve(&ct.masked_two_torsion_basis_EB.Q.x()),
        SUCCESS_RETVAL,
        "Q_B is not on the curve E_B",
    );
    assert_eq!(
        ct.codomain_curve
            .is_on_curve(&ct.masked_two_torsion_basis_EB.PQ.x()),
        SUCCESS_RETVAL,
        "P_B - Q_B is not on the curve E_B",
    );

    // Check (P_AB, Q_AB, P_AB - Q_AB), a basis of the 2^a-torsion on E_AB
    assert_eq!(
        ct.shared_end_curve
            .is_on_curve(&ct.masked_two_torsion_basis_EAB.P.x()),
        SUCCESS_RETVAL,
        "P_AB is not on the curve E_AB",
    );
    assert_eq!(
        ct.shared_end_curve
            .is_on_curve(&ct.masked_two_torsion_basis_EAB.Q.x()),
        SUCCESS_RETVAL,
        "Q_AB is not on the curve E_AB",
    );
    assert_eq!(
        ct.shared_end_curve
            .is_on_curve(&ct.masked_two_torsion_basis_EAB.PQ.x()),
        SUCCESS_RETVAL,
        "P_AB - Q_AB is not on the curve E_AB",
    );

    // Check (X_B, Y_B, X_B - Y_B), a basis of the 5^c-torsion on E_B
    assert_eq!(
        ct.codomain_curve
            .is_on_curve(&ct.masked_five_torsion_basis_EB.P.x()),
        SUCCESS_RETVAL,
        "X_B is not on the curve E_B",
    );
    assert_eq!(
        ct.codomain_curve
            .is_on_curve(&ct.masked_five_torsion_basis_EB.Q.x()),
        SUCCESS_RETVAL,
        "Y_B is not on the curve E_B",
    );
    assert_eq!(
        ct.codomain_curve
            .is_on_curve(&ct.masked_five_torsion_basis_EB.PQ.x()),
        SUCCESS_RETVAL,
        "X_B - Y_B is not on the curve E_B",
    );
}

#[apply(encryption_test_data)]
fn masked_image_points_are_linearly_independent<
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
    params: PublicParams<
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
    pub_key: PubKey<Fp2>,
) where
    [(); Fp2::ENCODED_LENGTH]: Sized,
{
    let mut message = *MESSAGE;
    let message = message.as_mut_slice();
    let (ct, _) = encrypt(&params, &pub_key, message);

    // Check (P_B, Q_B, P_B - Q_B), a basis of the 2^a-torsion on E_B
    let pair = ct.codomain_curve.weil_pairing_2exp(
        &ct.masked_two_torsion_basis_EB.P.x(),
        &ct.masked_two_torsion_basis_EB.Q.x(),
        &ct.masked_two_torsion_basis_EB.PQ.x(),
        params.two_torsion.exp,
    );
    assert_eq!(
        pair.pow(
            &params.two_torsion.reduced_order.to_le_bytes(),
            params.two_torsion.reduced_order.nbits(),
        )
        .equals(&Fp2::ONE),
        FAILURE_RETVAL,
    );
    assert_eq!(
        pair.pow(
            &params.two_torsion.order.to_le_bytes(),
            params.two_torsion.order.nbits(),
        )
        .equals(&Fp2::ONE),
        SUCCESS_RETVAL,
    );

    // Check (P_AB, Q_AB, P_AB - Q_AB), a basis of the 2^a-torsion on E_AB
    let pair = ct.shared_end_curve.weil_pairing(
        &ct.masked_two_torsion_basis_EAB.P.x(),
        &ct.masked_two_torsion_basis_EAB.Q.x(),
        &ct.masked_two_torsion_basis_EAB.PQ.x(),
        &params.two_torsion.order.to_le_bytes(),
        params.two_torsion.order.nbits(),
    );
    assert_eq!(
        pair.pow(
            &params.two_torsion.reduced_order.to_le_bytes(),
            params.two_torsion.reduced_order.nbits(),
        )
        .equals(&Fp2::ONE),
        FAILURE_RETVAL,
    );
    assert_eq!(
        pair.pow(
            &params.two_torsion.order.to_le_bytes(),
            params.two_torsion.order.nbits(),
        )
        .equals(&Fp2::ONE),
        SUCCESS_RETVAL,
    );

    // Check (X_B, Y_B, X_B - Y_B), a basis of the 5^c-torsion on E_B
    let pair = ct.codomain_curve.weil_pairing(
        &ct.masked_five_torsion_basis_EB.P.x(),
        &ct.masked_five_torsion_basis_EB.Q.x(),
        &ct.masked_five_torsion_basis_EB.PQ.x(),
        &params.five_torsion.order.to_le_bytes(),
        params.five_torsion.order.nbits(),
    );
    assert_eq!(
        pair.pow(
            &params.five_torsion.reduced_order.to_le_bytes(),
            params.five_torsion.reduced_order.nbits(),
        )
        .equals(&Fp2::ONE),
        FAILURE_RETVAL,
    );
    assert_eq!(
        pair.pow(
            &params.five_torsion.order.to_le_bytes(),
            params.five_torsion.order.nbits(),
        )
        .equals(&Fp2::ONE),
        SUCCESS_RETVAL,
    );
}
