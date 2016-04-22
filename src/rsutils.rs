
use std::mem;
    
pub fn double_to_raw_long_bits(double: &f64) -> i64 {
	unsafe { mem::transmute(*double) }	
}

pub fn long_bits_to_double(long: &i64) -> f64 {
	unsafe { mem::transmute(*long) }
}

pub fn float_to_raw_int_bits(float: &f32) -> i32 {
	unsafe { mem::transmute(*float) }	
}

pub fn int_bits_to_float(int: &i32) -> f32 {
	unsafe { mem::transmute(*int) }
}

#[cfg(test)]
mod tests {	
	use rsutils::*;
	
	#[test]
	fn double_to_long() {
		let double: f64 = -1.0;
		for a in 1..10000000 {
			let d = double + 1.0 / a as f64;
			assert!(d == long_bits_to_double(&double_to_raw_long_bits(&d)));
			let md = double - 1.0 / a as f64;
			assert!(md == long_bits_to_double(&double_to_raw_long_bits(&md)));
		}
	}

	#[test]
	fn float_to_long() {
		let float: f32 = -1.0;
		for a in 1..10000000 {
			let d = float + 1.0  / a as f32;
			assert!(d == int_bits_to_float(&float_to_raw_int_bits(&d)));
			let md = float - 1.0 / a as f32;
			assert!(md == int_bits_to_float(&float_to_raw_int_bits(&md)));
		}
	}
}