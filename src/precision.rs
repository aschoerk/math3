#![allow(dead_code)] 
#![allow(non_snake_case)]
use std::f64;
use rsutils::*;
use fastmath;

/* Exponent offset in IEEE754 representation. */
const EXPONENT_OFFSET: i64 = 1023; 
/* Offset to order signed double numbers lexicographically. */
const SGN_MASK: u64 = 0x8000000000000000;
/* Offset to order signed double numbers lexicographically. */
const SGN_MASK_FLOAT: u32 = 0x80000000;
/* Positive zero. */
const POSITIVE_ZERO: f64 = 0.0;

struct Constants {
	POSITIVE_ZERO_DOUBLE_BITS: i64,
	NEGATIVE_ZERO_DOUBLE_BITS: i64,
	POSITIVE_ZERO_FLOAT_BITS: i32,
	NEGATIVE_ZERO_FLOAT_BITS: i32,
	EPSILON: f64,
	SAFE_MIN: f64
}

impl Constants {
	fn new() -> Constants {
		Constants { 
			POSITIVE_ZERO_DOUBLE_BITS: double_to_raw_long_bits(&0.0),
			NEGATIVE_ZERO_DOUBLE_BITS: double_to_raw_long_bits(&-0.0),
			POSITIVE_ZERO_FLOAT_BITS:  float_to_raw_int_bits(&0.0),
			NEGATIVE_ZERO_FLOAT_BITS:  float_to_raw_int_bits(&-0.0),
			EPSILON: long_bits_to_double(&((EXPONENT_OFFSET - 53) << 52)),
			SAFE_MIN: long_bits_to_double(&((EXPONENT_OFFSET - 1022) << 52)),
		}
	}
	
	fn get() -> Constants {
		use std::sync:: {Once, ONCE_INIT};
		static START: Once = ONCE_INIT;
		static mut constants: *Constants = 0 as *Constants;
		unsafe {
			START.call_once(|| {
						constants = Constants::new();
			});
		};
	}
}

struct Double {
	
}

impl Double {
	
	fn equals_in_ulps(x: f64, y: f64, maxUlps: i32) {
		let xInt = double_to_raw_long_bits(x);
		let yInt = double_to_raw_long_bits(y);
		
		let isEqual =
		if ((xInt ^ yInt) & SGN_MASK) == 0 {
			fastmath::abs(xInt - yInt) <= maxUlps
		} else {
			let (deltaPlus, deltaMinus) =
			if xInt < yInt {
				(yInt - Constants::get().POSITIVE_ZERO_DOUBLE_BITS, xInt - Constants::get().NEGATIVE_ZERO_DOUBLE_BITS)
			} else {
				(xInt - Constants::get().POSITIVE_ZERO_DOUBLE_BITS, yInt - Constants::get().NEGATIVE_ZERO_DOUBLE_BITS)
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
		
	}

	fn equals_in_1(x: f64, y: f64) -> bool {
		true
	}
}

 