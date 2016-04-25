extern crate math3;

use math3::rsutils::*;

pub fn main() {
	
	
	let double: f64 = -1.0;
	for a in 1..10000000 {
		let d = double + 1.0 / a as f64;
		assert!(d == long_bits_to_double(&double_to_raw_long_bits(&d)));
		let md = double - 1.0 / a as f64;
		assert!(md == long_bits_to_double(&double_to_raw_long_bits(&md)));
	}

	
	println!("hello")
}