use poke::{
    SUCCESS_RETVAL,
    inke::{decrypt, encrypt, keygen},
    params,
};

fn main() {
    let params = params::inke_iii::get_params();

    let (pub_key, prv_key, ok) = keygen(&params);
    assert_eq!(ok, SUCCESS_RETVAL, "Key generation finished with errors");

    let message = String::from("Hello, world!");
    let message = message.as_bytes();
    println!("Original message: {:?}", message);

    let (ct, ok) = encrypt(&params, &pub_key, message);
    assert_eq!(ok, SUCCESS_RETVAL, "Encryption finished with errors");
    println!("Encrypted message: {:?}", ct.encrypted_message);

    let (message, ok) = decrypt(&params, &prv_key, &ct);
    assert_eq!(ok, SUCCESS_RETVAL, "Decryption finished with errors");
    println!("Decrypted message: {:?}", message);
}
