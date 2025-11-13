#![allow(non_snake_case)]

use poke::{SUCCESS_RETVAL, encrypt, example_keypairs, params};

fn main() {
    let params = params::poke_iii::get_params();
    let pub_key = example_keypairs::poke_iii::get_pub_key();

    let mut message = String::from("Hello, world!");
    let message = unsafe { message.as_bytes_mut() };
    let (ct, ok) = encrypt(&params, &pub_key, message);

    assert_eq!(ok, SUCCESS_RETVAL, "Encryption finished with errors");
    println!("Encrypted message: {:?}", ct.encrypted_message);
}
