use std::ops::Sub;
use std::cmp;
use std::cmp::Ordering;

pub struct Assert {
	
}

trait A {
	fn a(a: f64, b: f64) -> f64;
}

trait B {
	fn a(m: &str, a: f64, b: f64) -> f64;
}

impl Assert {
	
	fn differs_within_delta<T: PartialOrd + Sub<Output=T>>(expected: T, actual: T, delta: T) -> bool {
		Assert::differs_within_delta_m("differs_within_delta", expected, actual, delta)
	}
	
	fn differs_within_delta_m<T: PartialOrd + Sub<Output=T>>(message: &str, expected: T, actual: T, delta: T) -> bool {
		let cmp_res = expected.partial_cmp(&actual);
		match cmp_res {
			Some(Ordering::Equal) => true,
			Some(Ordering::Less) => delta > (actual - expected),
			Some(Ordering::Greater) => delta > (expected - actual),
			_ => false,
		}
	}
	
}

impl A for Assert {
	fn a(a: f64, b: f64) -> f64 {
		Assert::differs_within_delta(a,b,0.0000001);
		2.0 * a * b
	}
}

impl B for Assert {
	fn a(m: &str, a: f64, b: f64) -> f64 {
		Assert::differs_within_delta_m(m,a,b,0.0000001);
		1.0 * a
	}
}

#[cfg(test)]
mod tests {
	use assert;
	use assert::Assert;
	use assert::A;
	use assert::B;
	
	#[test]
	fn eqwd1() {
		assert!(Assert::differs_within_delta(1.1, 1.2, 0.2));
		assert!(<Assert as B>::a("delta",1.0,2.0) == 1.0);
		assert!(<Assert as A>::a(1.0,2.0) == 4.0);
	}
}