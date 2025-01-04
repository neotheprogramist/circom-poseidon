mod constants;
mod field;
mod poseidon;

use ark_ff::{Field, PrimeField};
use field::Fr;

fn main() {
    println!("Hello, world!");

    let x = Fr::from_be_bytes_mod_order(&prefix_hex::decode::<Vec<u8>>("0x01").unwrap());
    let y = Fr::from(2);

    let y = y.pow([5]);

    println!("x: {x} y: {y} x + y: {}", x + y);
}
