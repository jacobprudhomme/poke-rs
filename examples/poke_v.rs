use poke::{
    SUCCESS_RETVAL, example_keypairs, params,
    poke::{decrypt, encrypt},
};

fn main() {
    let params = params::poke_v::get_params();
    let pub_key = example_keypairs::poke_v::get_pub_key();
    let prv_key = example_keypairs::poke_v::get_prv_key();

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
