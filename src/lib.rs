#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]

pub mod bn;
mod dimtwo;
pub mod dlp;
pub mod example_keypairs;
pub mod fields;
pub mod inke;
mod masking;
mod modular;
pub mod params;
pub mod poke;
pub mod rand;

pub const SUCCESS_RETVAL: u32 = u32::MAX;
pub const FAILURE_RETVAL: u32 = u32::MIN;
