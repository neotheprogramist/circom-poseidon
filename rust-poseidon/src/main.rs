use field::Fr;

mod constants;
mod field;
mod functional;

fn main() {
    println!("Hello, world!");

    let result = functional::poseidon_hash(Fr::from(1), Fr::from(2));
    println!("result: {result}");
}
