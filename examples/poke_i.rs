#![allow(non_snake_case)]

use poke::{SUCCESS_RETVAL, decrypt, encrypt, example_keypairs, params};

fn main() {
    let params = params::poke_i::get_params();
    let pub_key = example_keypairs::poke_i::get_pub_key();
    let prv_key = example_keypairs::poke_i::get_prv_key();

    let mut message = String::from("Hello, world!");
    let message = unsafe { message.as_bytes_mut() };
    println!("Original message: {:?}", message);

    let (mut ct, ok) = encrypt(&params, &pub_key, message);
    assert_eq!(ok, SUCCESS_RETVAL, "Encryption finished with errors");
    println!("\nEncrypted message: {:?}", ct.encrypted_message);

    let (message, ok) = decrypt(&params, &prv_key, &mut ct);
    assert_eq!(ok, SUCCESS_RETVAL, "Decryption finished with errors");
    println!("\nDecrypted message: {:?}", message);
}
