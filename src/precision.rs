#![allow(dead_code)] 
#![allow(non_snake_case)]
#![allow(unused_variables)]
use std::f64;
use rsutils::*;
use fastmath;

/* Exponent offset in IEEE754 representation. */
const EXPONENT_OFFSET: u64 = 1023; 
/* Offset to order signed double numbers lexicographically. */
const SGN_MASK: u64 = 0x8000000000000000;
/* Offset to order signed double numbers lexicographically. */
const SGN_MASK_FLOAT: u32 = 0x80000000;
/* Positive zero. */
const POSITIVE_ZERO: f64 = 0.0;

lazy_static! {
	pub static ref POSITIVE_ZERO_DOUBLE_BITS: i64 = double_to_raw_long_bits(&0.0) as i64;
	pub static ref NEGATIVE_ZERO_DOUBLE_BITS: i64 = double_to_raw_long_bits(&-0.0) as i64;
	pub static ref POSITIVE_ZERO_FLOAT_BITS:  u32 = float_to_raw_int_bits(&0.0);
	pub static ref NEGATIVE_ZERO_FLOAT_BITS:  u32 = float_to_raw_int_bits(&-0.0);
	pub static ref EPSILON: f64 = long_bits_to_double(&((EXPONENT_OFFSET - 53) << 52));
	pub static ref SAFE_MIN: f64 = long_bits_to_double(&((EXPONENT_OFFSET - 1022) << 52));
}


struct Double {
	
}

impl Double {
	
	fn equals_in_ulps(x: f64, y: f64, maxUlps: i64) -> bool {
		let xInt = double_to_raw_long_bits(&x) as i64;
		let yInt = double_to_raw_long_bits(&y) as i64;
		
		let isEqual =
		if ((xInt as u64 ^ yInt as u64) & SGN_MASK) == 0 {
			fastmath::I64::abs(xInt - yInt) <= maxUlps
		} else {
			let (deltaPlus, deltaMinus) =
			if xInt < yInt {
				(yInt - *POSITIVE_ZERO_DOUBLE_BITS, xInt - *NEGATIVE_ZERO_DOUBLE_BITS)
			} else {
				(xInt - *POSITIVE_ZERO_DOUBLE_BITS, yInt - *NEGATIVE_ZERO_DOUBLE_BITS)
			};
			if deltaPlus > maxUlps {
				false
			} else {
				deltaMinus <= (maxUlps - deltaPlus)
			}
		};
		isEqual && !f64::is_nan(x) && !f64::is_nan(y)	
	}
	
	fn equals(x: f64, y: f64, eps: f64) -> bool {
		false
	}

	fn equals_in_1(x: f64, y: f64) -> bool {
		true
	}
}

 