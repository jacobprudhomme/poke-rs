#![allow(non_snake_case)]

use poke::{SUCCESS_RETVAL, decrypt, encrypt, example_keypairs, params::poke_i};

fn main() {
    let params = poke_i::get_params();
    let pub_key = example_keypairs::poke_i::get_pub_key();
    let prv_key = example_keypairs::poke_i::get_prv_key();

    let mut message = String::from("Hello, world!");
    let message = unsafe { message.as_bytes_mut() };
    let (ct, ok) = encrypt(&params, &pub_key, message);

    assert_eq!(ok, SUCCESS_RETVAL, "Encryption finished with errors");
    println!("Encrypted message: {:?}", ct.encrypted_message);
}
