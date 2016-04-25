#![allow(dead_code)] 
#[cfg(test)] extern crate hamcrest;

#[macro_use]
extern crate lazy_static;

#[cfg(test)] #[macro_use] pub mod assert;

pub mod base;
// pub mod complex;
pub mod rsutils;
pub mod fastmath;
// pub mod precision;


static a:i64 = (0x28be60db << 32) | 0x9391054a;

static b: f64 = 1.0 / 3.0;

pub fn test() -> i32 {
	if b == 1.0 {
		return 2;
	}
	if b == 2.0 {
		return 3;
	}
	4
}

#[cfg(test)] 
mod tests {
	
	static a:i64 = (0x28be60db << 32) | 0x9391054a;

	static b: f64 = 1.0 / 3.0;

	
	pub fn test() -> i32 {
		if b == 1.0 {
			return 2;
		}
		if b == 2.0 {
			return 3;
		}
		4
	}
	
	#[test]
	fn t2() {
		assert!(test() == 4);
	}
}