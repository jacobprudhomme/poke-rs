use poke::{BigNum, byte_bn_to_word_bn, word_bn_to_byte_bn};
use rand::RngCore as _;
use rstest::rstest;

#[rstest]
fn word_bn_is_inverse_of_byte_bn() {
    let mut rng = rand::rng();
    let mut bn_bytes = vec![0u8; 43];
    rng.fill_bytes(&mut bn_bytes);
    let last_byte = bn_bytes[42];
    let mut i = 0;
    let trailing_zero_bits = loop {
        if i >= 8 {
            break 0;
        }
        if last_byte & (1 << (8 - i - 1)) > 0 {
            break i;
        }
        i += 1;
    };
    let bn = BigNum {
        repr: bn_bytes,
        bitlen: 43 * 8,
    };

    let bn_words = byte_bn_to_word_bn(&bn);
    let bn_bytes_roundtrip = word_bn_to_byte_bn(&bn_words);

    assert_eq!(bn, bn_bytes_roundtrip);
}

#[rstest]
fn byte_bn_is_inverse_of_word_bn() {
    let mut rng = rand::rng();
    let mut bn_bytes = [0u8; 48];
    rng.fill_bytes(&mut bn_bytes);
    for i in 43..48 {
        bn_bytes[i] = 0;
    }
    let (_, bn_words, _) = unsafe { bn_bytes.align_to::<u64>() };

    let bn = word_bn_to_byte_bn(&bn_words);
    let bn_words_roundtrip = byte_bn_to_word_bn(&bn);

    assert_eq!(bn_words, bn_words_roundtrip);
}
