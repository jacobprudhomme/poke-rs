use poke::{byte_bn_to_word_bn, word_bn_to_byte_bn};
use rand::RngCore as _;
use rstest::rstest;

#[rstest]
fn word_bn_is_inverse_of_byte_bn() {
    let mut rng = rand::rng();
    let mut num_bytes = [0u8; 43];
    rng.fill_bytes(&mut num_bytes);

    let num_words = byte_bn_to_word_bn(&num_bytes);
    let num_bytes_roundtrip = word_bn_to_byte_bn(&num_words);

    assert_eq!(num_bytes.as_slice(), num_bytes_roundtrip);
}

#[rstest]
fn byte_bn_is_inverse_of_word_bn() {
    let mut rng = rand::rng();
    let mut num_bytes = [0u8; 48];
    rng.fill_bytes(&mut num_bytes);
    for i in 43..48 {
        num_bytes[i] = 0;
    }
    let (_, num_words, _) = unsafe { num_bytes.align_to::<u64>() };

    let num_bytes = word_bn_to_byte_bn(&num_words);
    let num_words_roundtrip = byte_bn_to_word_bn(&num_bytes);

    assert_eq!(num_words, num_words_roundtrip);
}
