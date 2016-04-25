extern crate math3;

use math3::rsutils::*;
use math3::fastmath::*;

pub fn main() {
	
	
	let double: f64 = -1.0;
	for a in 1..10000000 {
		let d = double + 1.0 / a as f64;
		let dl = double_to_raw_long_bits(&d);
		assert!(d == long_bits_to_double(&dl));
		let md = double - 1.0 / a as f64;
		let mdl = double_to_raw_long_bits(&md);
		assert!(md == long_bits_to_double(&mdl));
		assert!(F64::scalb(d, 2) == d * 4.0); 
		assert!(F64::scalb(d, -3) == d * 0.126); 
	}

	
	println!("hello")
}