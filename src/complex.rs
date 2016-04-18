

use std::f64;

pub struct Complex {
	real: f64,
	
	imaginary: f64,
	
	is_nan: bool,
	
	is_infinite: bool,
	
}

pub const I: Complex = Complex {real: 0.0, imaginary: 1.0, is_nan: false, is_infinite: false };

    /** A complex number representing "NaN + NaNi" */
pub const NAN: Complex = Complex {real: f64::NAN, imaginary: f64::NAN, is_nan: true, is_infinite: false };

    /** A complex number representing "+INF + INFi" */
pub const INF: Complex = Complex {real: f64::INFINITY, imaginary: f64::INFINITY, is_nan: false, is_infinite: true };

    /** A complex number representing "1.0 + 0.0i" */
pub const ONE: Complex = Complex {real: 1.0, imaginary: 0.0, is_nan: false, is_infinite: false };
    
    /** A complex number representing "0.0 + 0.0i" */    
pub const ZERO: Complex = Complex {real: 0.0, imaginary: 0.0, is_nan: false, is_infinite: false };

impl PartialEq for Complex {
	 fn eq(&self, other: &Complex) -> bool {
	 	if other.is_nan {
	 		self.is_nan
	 	} else {
	 		self.real == other.real && self.imaginary == other.imaginary
	 	}
	 }
}

impl Complex {

	
	pub fn new(real: f64, imaginary: f64) -> Complex {
		let is_nan = f64::is_nan(real) || f64::is_nan(imaginary);
		Complex { real: real, imaginary: imaginary, 
			is_nan: is_nan, is_infinite: f64::is_infinite(real) || f64::is_infinite(imaginary)}
	}

	pub fn new_real(real: f64) -> Complex {
		Complex::new(real, 0.0)
	}
	
	pub fn get_real(&self) -> f64 {
		self.real
	}
	
	pub fn get_imaginary(&self) -> f64 {
		self.imaginary
	}
	
	pub fn is_nan(&self) -> bool {
		self.is_nan
	}
	
	pub fn is_infinite(&self) -> bool {
		self.is_infinite
	}
	
	
	pub fn abs(&self) -> f64 {
		match self {
			&Complex {real: _, imaginary:_,is_nan: true, is_infinite: _} => f64::NAN,
			&Complex {real: _, imaginary:_,is_nan: _, is_infinite: true} => f64::INFINITY,
			_ => if f64::abs(self.real) < f64::abs(self.imaginary) {
				if self.imaginary == 0.0 {
					f64::abs(self.real)
				} else {
					let q = self.real / self.imaginary;
					f64::abs(self.imaginary) * f64::sqrt(1.0 + q * q)
				}
			} else {
				if self.real == 0.0 {
					f64::abs(self.imaginary)
				} else {
					let q = self.imaginary / self.real;
					f64::abs(self.real) * f64::sqrt(1.0 + q * q)
				}
			},
		}
	}
	
	pub fn add(&self, addend: &Complex) -> Complex  {
        if self.is_nan || addend.is_nan {
            NAN
        } else {
          Complex::new (self.real + addend.get_real(),
                             self.imaginary + addend.get_imaginary())
        }
    }
	
	 pub fn add_f64(&self, addend: f64) -> Complex {
        if self.is_nan || f64::is_nan(addend) {
            NAN
        } else {
        	Complex::new(self.real + addend, self.imaginary)
        }
    }
	 
	 pub fn conjugate(&self) -> Complex {
        if self.is_nan {
            NAN
        } else {
          	Complex::new(self.real, -self.imaginary)
        }
    }
	 
	 pub fn divide(&self, divisor: &Complex) -> Complex {
        if self.is_nan || divisor.is_nan {
            return NAN
        }

        let c = divisor.get_real();
        let d = divisor.get_imaginary();
        if c == 0.0 && d == 0.0 {
            return NAN
        }

        if divisor.is_infinite() && !self.is_infinite() {
            return ZERO
        }

        if f64::abs(c) < f64::abs(d) {
            let q = c / d;
            let denominator = c * q + d;
            Complex::new((self.real * q + self.imaginary) / denominator,
                (self.imaginary * q - self.real) / denominator)
        } else {
            let q = d / c;
            let denominator = d * q + c;
            Complex::new((self.imaginary * q + self.real) / denominator,
                (self.imaginary - self.real * q) / denominator)
        }
    }
	 
	 pub fn divide_f64(&self, divisor: f64) -> Complex {
        if self.is_nan || f64::is_nan(divisor) {
            return NAN
        }
        if divisor == 0.0 {
            return NAN
        }
        if f64::is_infinite(divisor) {
            if !self.is_infinite { ZERO } else { NAN }
        } else {
        	Complex::new(self.real / divisor,
                             self.imaginary  / divisor)
        }
    }
	 
	 pub fn reciprocal(&self) -> Complex {
        if self.is_nan {
            return NAN
        }

        if self.real == 0.0 && self.imaginary == 0.0 {
            return INF
        }

        if self.is_infinite {
            return ZERO
        }

        if f64::abs(self.real) < f64::abs(self.imaginary) {
            let q = self.real / self.imaginary;
            let scale = 1.0 / (self.real * q + self.imaginary);
            Complex::new(scale * q, -scale)
        } else {
            let q = self.imaginary / self.real;
            let scale = 1. / (self.imaginary * q + self.real);
            Complex::new(scale, -scale * q)
        }
    }
}


#[cfg(test)]
mod tests {
	use hamcrest::*;	
	use complex;
	use complex::Complex;	
	use std::f64;
	
	const ONE_INF: Complex = Complex {real: 1.0, imaginary: f64::INFINITY, is_nan: false, is_infinite: true};
	const ONE_NEG_INF: Complex = Complex {real: 1.0, imaginary: f64::NEG_INFINITY, is_nan: false, is_infinite: true};
	const INF_ONE: Complex = Complex {real: f64::INFINITY, imaginary: 1.0, is_nan: false, is_infinite: true};
	const INF_ZERO: Complex = Complex {real: f64::INFINITY , imaginary: 0.0, is_nan: false, is_infinite: true};
	const INF_NAN: Complex = Complex {real: f64::INFINITY, imaginary: f64::NAN, is_nan: true, is_infinite: true};
	const INF_NEG_INF: Complex = Complex {real: f64::INFINITY, imaginary: f64::NEG_INFINITY, is_nan: false, is_infinite: true};
	const NEG_INF_INF: Complex = Complex {real: f64::NEG_INFINITY, imaginary: f64::INFINITY, is_nan: false, is_infinite: true};
	const NEG_INF_ZERO: Complex = Complex {real: f64::NEG_INFINITY, imaginary: 0.0, is_nan: false, is_infinite: true};
	const NEG_INF_ONE: Complex = Complex {real: f64::NEG_INFINITY, imaginary: 1.0, is_nan: false, is_infinite: true};
	const NEG_INF_NAN: Complex = Complex {real: f64::NEG_INFINITY, imaginary: f64::NAN, is_nan: true, is_infinite: true};
	const NEG_INF_NEG_INF: Complex = Complex {real: f64::NEG_INFINITY, imaginary: f64::NEG_INFINITY, is_nan: false, is_infinite: true};
	const ONE_NAN: Complex = Complex {real: 1.0, imaginary: f64::NAN, is_nan: true, is_infinite: false};
	const ZERO_INF: Complex = Complex {real: 0.0, imaginary: f64::INFINITY, is_nan: false, is_infinite: true};
	const ZERO_NAN: Complex = Complex {real: 0.0, imaginary: f64::NAN, is_nan: true, is_infinite: false};
	const NAN_INF: Complex = Complex {real: f64::NAN, imaginary: f64::INFINITY, is_nan: true, is_infinite: true};
	const NAN_NEG_INF: Complex = Complex {real: f64::NAN, imaginary: f64::NEG_INFINITY, is_nan: true, is_infinite: true};
	const NAN_ZERO: Complex = Complex {real: f64::NAN, imaginary: 0.0, is_nan: true, is_infinite: false};
    
	#[test]
	fn nan() {
		let c: Complex = Complex::new(f64::NAN, f64::NAN);
		assert!(c.is_nan);
	}
	
	#[test]
	fn constructor_test() {
		let c: Complex = Complex::new(3.0, 4.0);
		assert_that(c.real,is(close_to(3.0,0.00001)));
		assert_that(c.imaginary,is(close_to(4.0,0.00001)));
		
	}
	
	#[test]
	fn test_constructor_nan() {
        let z: Complex = Complex::new(3.0, f64::NAN);
        assert!(z.is_nan);

        let z: Complex = Complex::new(f64::NAN, 4.0);
        assert!(z.is_nan);
 
        let z: Complex = Complex::new(3.0, 4.0);
        assert!(!z.is_nan);
    }
	
	#[test]
	fn test_abs() {
		let z: Complex = Complex::new(3.0,4.0);
		assert_that(z.abs(), is(close_to(5.0,0.0001)))
	}
	
	#[test]
    fn test_abs_nan() {
        assert!(f64::is_nan(complex::NAN.abs()));
        let z: Complex = Complex::new(f64::INFINITY, f64::NAN);
        assert!(f64::is_nan(z.abs()));
    }

    #[test]
    fn test_abs_infinite() {
        let z: Complex = Complex::new(f64::INFINITY, 0.0);
        assert_that(f64::INFINITY,is(close_to( z.abs(), 0.0)));
        let z: Complex = Complex::new(0.0, f64::NEG_INFINITY);
        assert_that(f64::INFINITY,is(close_to( z.abs(), 0.0)));
        let z: Complex = Complex::new(f64::INFINITY, f64::NEG_INFINITY);
        assert_that(f64::INFINITY,is(close_to( z.abs(), 0.0)));
    }

    #[test]
    fn test_add() {
        let x: Complex = Complex::new(3.0, 4.0);
        let y: Complex = Complex::new(5.0, 6.0);
        let z: Complex = x.add(&y);
        assert_that(8.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(10.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testadd_nan() {
        let x: Complex = Complex::new(3.0, 4.0);
        let z: Complex = x.add(&complex::NAN);
        assert!(z.is_nan());
        let z = Complex::new(1.0, f64::NAN);
        let w: Complex = x.add(&z);
        assert!(w.is_nan());
    }

    #[test]
    fn test_add_inf() {
        let x = Complex::new(1.0, 1.0);
        let z = Complex::new(f64::INFINITY, 0.0);
        let w: Complex = x.add(&z);
        assert!(w.get_imaginary() ==  1.0);
        assert!(f64::INFINITY == w.get_real());

        let x = Complex::new(f64::NEG_INFINITY, 0.0);
        assert!(f64::is_nan(x.add(&z).get_real()));
    }


    #[test]
    fn test_scalar_add() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = 2.0;
        let y_complex = Complex::new_real(y_double);
        assert!(x.add(&y_complex) == x.add_f64(y_double));
    }

    #[test]
    fn test_scalar_add_nan() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = f64::NAN;
        let y_complex = Complex::new_real(y_double);
        assert!(x.add(&y_complex) == x.add_f64(y_double));
    }

    #[test]
    fn test_scalar_add_inf() {
        let x = Complex::new(1.0, 1.0);
        let y_double: f64 = f64::INFINITY;

        let y_complex = Complex::new_real(y_double);
        assert!(x.add(&y_complex) == x.add_f64(y_double));

        let x = Complex::new(f64::INFINITY, 0.0);
        assert!(x.add(&y_complex) == x.add_f64(y_double));
    }


    #[test]
    fn test_conjugate() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.conjugate();
        assert_that(3.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(-4.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn test_conjugate_nan() {
        let z = complex::NAN.conjugate();
        assert!(z.is_nan());
    }


    #[test]
    fn test_conjugate_infiinite() {
        let z = Complex::new(0.0, f64::INFINITY);
        assert_that(f64::INFINITY, is(close_to(z.conjugate().get_imaginary(), 0.0)));
        let z = Complex::new(0.0, f64::INFINITY);
        assert_that(f64::INFINITY, is(close_to(z.conjugate().get_imaginary(), 0.0)));
    }

    #[test]
    fn test_divide() {
        let x = Complex::new(3.0, 4.0);
        let y = Complex::new(5.0, 6.0);
        let z: Complex = x.divide(&y);
        assert_that(39.0 / 61.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(2.0 / 61.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn test_divide_real() {
        let x = Complex::new(2f64, 3f64);
        let y = Complex::new(2f64, 0f64);
        assert!(Complex::new(1f64, 1.5f64) == x.divide(&y));

    }


    #[test]
    fn test_divide_imaginary() {
        let x = Complex::new(2f64, 3f64);
        let y = Complex::new(0f64, 2f64);
        assert!(Complex::new(1.5f64	, -1f64) == x.divide(&y));
    }

    #[test]
    fn test_divide_inf() {
        let x = Complex::new(3.0, 4.0);
        let w = Complex::new(f64::NEG_INFINITY, f64::INFINITY);
        assert!(x.divide(&w) == complex::ZERO);

        let z = w.divide(&x);
        assert!(f64::is_nan(z.get_real()));
        assert_that(f64::INFINITY, is(close_to(z.get_imaginary(), 0.0)));

        let w = Complex::new(f64::INFINITY, f64::INFINITY);
        let z = w.divide(&x);
        assert!(f64::is_nan(z.get_imaginary()));
        assert_that(f64::INFINITY, is(close_to(z.get_real(), 0.0)));

        let w = Complex::new(1.0, f64::INFINITY);
        let z = w.divide(&w);
        assert!(f64::is_nan(z.get_real()));
        assert!(f64::is_nan(z.get_imaginary()));
    }

    #[test]
    fn test_divide_zero() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.divide(&complex::ZERO);
        // assert_that(z, Complex.INF); // See MATH-657
        assert!(z == complex::NAN);
    }

    #[test]
    fn test_divide_zero_zero() {
        let x = Complex::new(0.0, 0.0);
        let z: Complex = x.divide(&complex::ZERO);
        assert!(z == complex::NAN);
    }

    #[test]
    fn test_divide_nan() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.divide(&complex::NAN);
        assert!(z.is_nan());
    }

    #[test]
    fn test_divide_nan_inf() {
       let z: Complex = ONE_INF.divide(&complex::ONE);
       assert!(f64::is_nan(z.get_real()));
       assert_that(f64::INFINITY,is(close_to( z.get_imaginary(), 0.0)));

       let z = NEG_INF_NEG_INF.divide(&ONE_NAN);
       assert!(f64::is_nan(z.get_real()));
       assert!(f64::is_nan(z.get_imaginary()));

       let z = NEG_INF_INF.divide(&complex::ONE);
       assert!(f64::is_nan(z.get_real()));
       assert!(f64::is_nan(z.get_imaginary()));
    }

    #[test]
    fn test_scalar_divide() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = 2.0;
        let y_complex = Complex::new_real(y_double);
        assert!(x.divide(&y_complex) == x.divide_f64(y_double));
    }

    #[test]
    fn test_scalar_divide_nan() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = f64::NAN;
        let y_complex = Complex::new_real(y_double);
        assert!(x.divide(&y_complex) == x.divide_f64(y_double));
    }

    #[test]
    fn test_scalar_divide_inf() {
        let x = Complex::new(1.0,1.0);
        let y_double: f64 = f64::INFINITY;
        let y_complex = Complex::new_real(y_double);
        assert!(x.divide(&y_complex) == x.divide_f64(y_double));

        let y_double = f64::NEG_INFINITY;
        let y_complex = Complex::new_real(y_double);
        assert!(x.divide(&y_complex) == x.divide_f64(y_double));

        let x = Complex::new(1.0, f64::NEG_INFINITY);
        assert!(x.divide(&y_complex) == x.divide_f64(y_double));
    }

    #[test]
    fn test_scalarDivideZero() {
        let x = Complex::new(1.0,1.0);
        assert!(x.divide(complex::ZERO) == x.divide_f64(0.0));
    }

    #[test]
    fn testReciprocal() {
        let z = Complex::new(5.0, 6.0);
        let act = z.reciprocal();
        let expRe = 5.0 / 61.0;
        let expIm = -6.0 / 61.0;
        assert_that(expRe, is(close_to(act.get_real(), FastMath.ulp(expRe))));
        assert_that(expIm, is(close_to(act.get_imaginary(), FastMath.ulp(expIm))));
    }

/*
    #[test]
    fn testReciprocalReal() {
        let z = Complex::new(-2.0, 0.0);
        assert!(Complex.equals(Complex::new(-0.5, 0.0), z.reciprocal()));
    }

    #[test]
    fn testReciprocalImaginary() {
        let z = Complex::new(0.0, -2.0);
        assert_that(Complex::new(0.0, 0.5), z.reciprocal());
    }

    #[test]
    fn testReciprocalInf() {
        let z = Complex::new(f64::INFINITY, f64::INFINITY);
        assert!(z.reciprocal().equals(complex::ZERO));

        z = Complex::new(1, f64::INFINITY).reciprocal();
        assert_that(z, complex::ZERO);
    }

    #[test]
    fn testReciprocalZero() {
        assert_that(complex::ZERO.reciprocal(), Complex.INF);
    }

    #[test]
    fn testReciprocal_nan() {
        assert!(complex::NAN.reciprocal().is_nan());
    }

    #[test]
    fn testMultiply() {
        let x = Complex::new(3.0, 4.0);
        let y = Complex::new(5.0, 6.0);
        let z: Complex = x.multiply(y);
        assert_that(-9.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(38.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testMultiply_nan() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.multiply(complex::NAN);
        Assert.assertSame(complex::NAN, z);
        z = complex::NAN.multiply(5);
        Assert.assertSame(complex::NAN, z);
    }

    #[test]
    fn testMultiplyInfInf() {
        // assert!(infInf.multiply(infInf).is_nan()); // MATH-620
        assert!(infInf.multiply(infInf).isInfinite());
    }

    #[test]
    fn testMultiply_nanInf() {
        let z = Complex::new(1,1);
        let w: Complex = z.multiply(infOne);
        assert_that(w.get_real(), is(close_to(f64::INFINITY, 0)));
        assert_that(w.get_imaginary(), is(close_to(f64::INFINITY, 0)));

        // [MATH-164]
        assert!(Complex::new( 1,0).multiply(infInf).equals(Complex.INF));
        assert!(Complex::new(-1,0).multiply(infInf).equals(Complex.INF));
        assert!(Complex::new( 1,0).multiply(negInfZero).equals(Complex.INF));

        w = oneInf.multiply(oneNegInf);
        assert_that(w.get_real(), is(close_to(f64::INFINITY, 0)));
        assert_that(w.get_imaginary(), is(close_to(f64::INFINITY, 0)));

        w = negInfNegInf.multiply(one_nan);
        assert!(f64::is_nan(w.get_real()));
        assert!(f64::is_nan(w.get_imaginary()));

        z = Complex::new(1, f64::INFINITY);
        Assert.assertSame(Complex.INF, z.multiply(z));
    }

    #[test]
    fn test_scalarMultiply() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = 2.0;
        let y_complex = Complex::new(y_double);
        assert_that(x.multiply(y_complex), x.multiply(y_double));
        zInt = -5;
        let z_complex = Complex::new(zInt);
        assert_that(x.multiply(z_complex), x.multiply(zInt));
    }

    #[test]
    fn test_scalarMultiply_nan() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = f64::NAN;
        let y_complex = Complex::new(y_double);
        assert_that(x.multiply(y_complex), x.multiply(y_double));
    }

    #[test]
    fn test_scalarMultiplyInf() {
        let x = Complex::new(1, 1);
        let y_double: f64 = f64::INFINITY;
        let y_complex = Complex::new(y_double);
        assert_that(x.multiply(y_complex), x.multiply(y_double));

        y_double = Double.NEGATIVE_INFINITY;
        y_complex = Complex::new(y_double);
        assert_that(x.multiply(y_complex), x.multiply(y_double));
    }

    #[test]
    fn testNegate() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.negate();
        assert_that(-3.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(-4.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testNegate_nan() {
        let z: Complex = complex::NAN.negate();
        assert!(z.is_nan());
    }

    #[test]
    fn test_subtract() {
        let x = Complex::new(3.0, 4.0);
        let y = Complex::new(5.0, 6.0);
        let z: Complex = x.subtract(y);
        assert_that(-2.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(-2.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn test_subtract_nan() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.subtract(complex::NAN);
        Assert.assertSame(complex::NAN, z);
        z = Complex::new(1, f64::NAN);
        let w: Complex = x.subtract(z);
        Assert.assertSame(complex::NAN, w);
    }

    #[test]
    fn test_subtractInf() {
        let x = Complex::new(1, 1);
        let z = Complex::new(f64::INFINITY, 0);
        let w: Complex = x.subtract(z);
        assert_that(w.get_imaginary(), is(close_to(1, 0)));
        assert_that(f64::INFINITY, is(close_to(w.get_real(), 0)));

        x = Complex::new(f64::INFINITY, 0);
        assert!(f64::is_nan(x.subtract(z).get_real()));
    }

    #[test]
    fn test_scalarSubtract() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = 2.0;
        let y_complex = Complex::new(y_double);
        assert_that(x.subtract(y_complex), x.subtract(y_double));
    }

    #[test]
    fn test_scalarSubtract_nan() {
        let x = Complex::new(3.0, 4.0);
        let y_double: f64 = f64::NAN;
        let y_complex = Complex::new(y_double);
        assert_that(x.subtract(y_complex), x.subtract(y_double));
    }

    #[test]
    fn test_scalarSubtractInf() {
        let x = Complex::new(1, 1);
        let y_double: f64 = f64::INFINITY;
        let y_complex = Complex::new(y_double);
        assert_that(x.subtract(y_complex), x.subtract(y_double));

        x = Complex::new(f64::INFINITY, 0);
        assert_that(x.subtract(y_complex), x.subtract(y_double));
    }
    */
}


