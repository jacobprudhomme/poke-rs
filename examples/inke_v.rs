use poke::{
    SUCCESS_RETVAL,
    inke::{decrypt, encrypt, keygen},
    params,
};
use rand::RngCore as _;

fn main() {
    let params = params::inke_v::get_params();

    let (pub_key, prv_key, ok) = keygen(&params);
    assert_eq!(ok, SUCCESS_RETVAL, "Key generation finished with errors");

    let mut rng = rand::rng();
    let mut message = [0; 32];
    rng.fill_bytes(&mut message);
    println!("Original message:\n{:?}\n", message);

    let (ct, ok) = encrypt(&params, &pub_key, &message);
    assert_eq!(ok, SUCCESS_RETVAL, "Encryption finished with errors");
    println!("Encrypted message:\n{:?}\n", ct.encrypted_message);

    let (dec_message, ok) = decrypt(&params, &prv_key, &ct);
    assert_eq!(ok, SUCCESS_RETVAL, "Decryption finished with errors");
    println!("Decrypted message:\n{:?}\n", dec_message);

    if dec_message == message {
        println!("SUCCESS: Messages match!");
    } else {
        println!("FAILURE: Messages do not match.");
    }
}
