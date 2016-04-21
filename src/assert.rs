use std::ops::Sub;
use std::cmp;
use std::cmp::Ordering;


pub struct Assert {
	
}

impl Assert {
		
	pub fn differs_within_delta<T: PartialOrd + Sub<Output=T>>(expected: T, actual: T, delta: T) -> bool {
		let cmp_res = expected.partial_cmp(&actual);
		match cmp_res {
			Some(Ordering::Equal) => true,
			Some(Ordering::Less) => delta > (actual - expected),
			Some(Ordering::Greater) => delta > (expected - actual),
			_ => false,
		}
	}
	
}

#[macro_export]
macro_rules! assert_within_delta{
	($a:expr, $b:expr, $d:expr) => ({
			let (a, b, d) = (&$a, &$b, &$d);
			assert!(assert::Assert::differs_within_delta(*a, *b, *d),
				"{} differs from {} by more than {}", *a, *b, *d);			
	})
}

#[cfg(test)]
mod tests {
	#[macro_use]
	use assert;
	use assert::Assert;
	
	#[test]
	fn eqwd1() {
		assert!(Assert::differs_within_delta(1.1, 1.2, 0.2));
		assert_within_delta!(1.1, 1.2, 0.2);
	}
}