

use std::f64;

pub struct Complex {
	pub real: f64,
	
	pub imaginary: f64,
	
	pub is_nan: bool,
	
	pub is_infinite: bool,
	
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


impl Complex {

	
	pub fn new(real: f64, imaginary: f64) -> Complex {
		let is_nan = f64::is_nan(real) || f64::is_nan(imaginary);
		Complex { real: real, imaginary: imaginary, 
			is_nan: is_nan, is_infinite: f64::is_infinite(real) || f64::is_infinite(imaginary)}
	}
	
	fn get_real(&self) -> f64 {
		self.real
	}
	
	fn get_imaginary(&self) -> f64 {
		self.imaginary
	}
	
	fn is_nan(&self) -> bool {
		self.is_nan
	}
	
	fn is_infinite(&self) -> bool {
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
	
	fn add(&self, addend: Complex) -> Complex  {
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
	 
	 pub fn divide(&self, divisor: Complex) -> Complex {
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
	use complex::Complex;	
	use std::f64;
	#[test]
	fn nan() {
		let c = Complex::new(f64::NAN, f64::NAN);
		assert!(c.is_nan);
	}
	
	#[test]
	fn constructor_test() {
		let c = Complex::new(3.0, 4.0);
		assert_that(c.real,is(close_to(3.0,0.00001)));
		assert_that(c.imaginary,is(close_to(4.0,0.00001)));
		
	}
	
	#[test]
	fn test_constructor_nan() {
        let z = Complex::new(3.0, f64::NAN);
        assert!(z.is_nan);

        let z = Complex::new(f64::NAN, 4.0);
        assert!(z.is_nan);
 
        let z = Complex::new(3.0, 4.0);
        assert!(!z.is_nan);
    }
	
	#[test]
	fn test_abs() {
		let z = Complex::new(3.0,4.0);
		assert_that(z.abs(), is(close_to(5.0,0.0001)))
	}
	
	#[test]
    fn testAbsNaN() {
        assert!(f64::is_nan(Complex::NAN.abs()));
        let z = Complex::new(f64::INFINITY, f64::NAN);
        assert!(f64::is_nan(z.abs()));
    }

    #[test]
    fn testAbsInfinite() {
        let z = Complex::new(f64::INFINITY, 0);
        assert_that(f64::INFINITY,is(close_to( z.abs(), 0)));
        z = Complex::new(0, f64::INFINITY);
        assert_that(f64::INFINITY,is(close_to( z.abs(), 0)));
        z = Complex::new(f64::INFINITY, f64::INFINITY);
        assert_that(f64::INFINITY,is(close_to( z.abs(), 0)));
    }

    #[test]
    fn testAdd() {
        let x = Complex::new(3.0, 4.0);
        let y = Complex::new(5.0, 6.0);
        let z: Complex = x.add(y);
        assert_that(8.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(10.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testAddNaN() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.add(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = Complex::new(1, f64::NAN);
        let w: Complex = x.add(z);
        Assert.assertSame(Complex.NaN, w);
    }

    #[test]
    fn testAddInf() {
        let x = Complex::new(1, 1);
        let z = Complex::new(f64::INFINITY, 0);
        let w: Complex = x.add(z);
        assert_that(w.get_imaginary(),is(close_to( 1, 0)));
        assert_that(f64::INFINITY, is(close_to(w.get_real(), 0)));

        x = Complex::new(f64::INFINITY, 0);
        assert!(f64::is_nan(x.add(z).get_real()));
    }


    #[test]
    fn testScalarAdd() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = 2.0;
        let yComplex = Complex::new(yDouble);
        assert_that(x.add(yComplex),is(close_to( x.add(yDouble))));
    }

    #[test]
    fn testScalarAddNaN() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = Double.NaN;
        let yComplex = Complex::new(yDouble);
        assert_that(x.add(yComplex),is(close_to( x.add(yDouble))));
    }

    #[test]
    fn testScalarAddInf() {
        let x = Complex::new(1, 1);
        let yDouble: f64 = Double.POSITIVE_INFINITY;

        let yComplex = Complex::new(yDouble);
        assert_that(x.add(yComplex),is(close_to( x.add(yDouble))));

        x = Complex::new(f64::INFINITY, 0);
        assert_that(x.add(yComplex),is(close_to( x.add(yDouble))));
    }

    #[test]
    fn testConjugate() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.conjugate();
        assert_that(3.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(-4.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testConjugateNaN() {
        let z = Complex.NaN.conjugate();
        assert!(z.isNaN());
    }

    #[test]
    fn testConjugateInfiinite() {
        let z = Complex::new(0, f64::INFINITY);
        assert_that(f64::INFINITY, is(close_to(z.conjugate().get_imaginary(), 0)));
        z = Complex::new(0, f64::INFINITY);
        assert_that(f64::INFINITY, is(close_to(z.conjugate().get_imaginary(), 0)));
    }

    #[test]
    fn testDivide() {
        let x = Complex::new(3.0, 4.0);
        let y = Complex::new(5.0, 6.0);
        let z: Complex = x.divide(y);
        assert_that(39.0 / 61.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(2.0 / 61.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testDivideReal() {
        let x = Complex::new(2f64, 3f64);
        let y = Complex::new(2f64, 0f64);
        assert!(Complex::new(1f64, 1.5f64) == x.divide(y));

    }

    #[test]
    fn testDivideImaginary() {
        let x = Complex::new(2f64, 3f64);
        let y = Complex::new(0f64, 2f64);
        assert!(Complex::new(1.5f64	, -1f64) == x.divide(y));
    }

    #[test]
    fn testDivideInf() {
        let x = Complex::new(3, 4);
        let w = Complex::new(f64::INFINITY, f64::INFINITY);
        assert!(x.divide(w).equals(Complex.ZERO));

        let z = w.divide(x);
        assert!(f64::is_nan(z.get_real()));
        assert_that(f64::INFINITY, z.get_imaginary(), 0);

        let w = Complex::new(f64::INFINITY, f64::INFINITY);
        let z = w.divide(x);
        assert!(f64::is_nan(z.get_imaginary()));
        assert_that(f64::INFINITY, z.get_real(), 0);

        let w = Complex::new(1, f64::INFINITY);
        let z = w.divide(w);
        assert!(f64::is_nan(z.get_real()));
        assert!(f64::is_nan(z.get_imaginary()));
    }

    #[test]
    fn testDivideZero() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.divide(Complex.ZERO);
        // assert_that(z, Complex.INF); // See MATH-657
        assert!(z == Complex.NaN);
    }

    #[test]
    fn testDivideZeroZero() {
        let x = Complex::new(0.0, 0.0);
        let z: Complex = x.divide(Complex.ZERO);
        assert!(z == Complex.NaN);
    }

    #[test]
    fn testDivideNaN() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.divide(Complex.NaN);
        assert!(z.isNaN());
    }

    #[test]
    fn testDivideNaNInf() {
       let z: Complex = oneInf.divide(Complex.ONE);
       assert!(f64::is_nan(z.get_real()));
       assert_that(f64::INFINITY,is(close_to( z.get_imaginary(), 0)));

       z = negInfNegInf.divide(oneNaN);
       assert!(f64::is_nan(z.get_real()));
       assert!(f64::is_nan(z.get_imaginary()));

       z = negInfInf.divide(Complex.ONE);
       assert!(f64::is_nan(z.get_real()));
       assert!(f64::is_nan(z.get_imaginary()));
    }

    #[test]
    fn testScalarDivide() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = 2.0;
        let yComplex = Complex::new(yDouble);
        assert!(x.divide(yComplex) == x.divide(yDouble));
    }

    #[test]
    fn testScalarDivideNaN() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = Double.NaN;
        let yComplex = Complex::new(yDouble);
        assert!(x.divide(yComplex) == x.divide(yDouble));
    }

    #[test]
    fn testScalarDivideInf() {
        let x = Complex::new(1,1);
        let yDouble: f64 = Double.POSITIVE_INFINITY;
        let yComplex = Complex::new(yDouble);
        TestUtils.assertEquals(x.divide(yComplex), x.divide(yDouble), 0);

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = Complex::new(yDouble);
        TestUtils.assertEquals(x.divide(yComplex), x.divide(yDouble), 0);

        x = Complex::new(1, Double.NEGATIVE_INFINITY);
        TestUtils.assertEquals(x.divide(yComplex), x.divide(yDouble), 0);
    }

    #[test]
    fn testScalarDivideZero() {
        let x = Complex::new(1,1);
        TestUtils.assertEquals(x.divide(Complex.ZERO), x.divide(0), 0);
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
        assert!(z.reciprocal().equals(Complex.ZERO));

        z = Complex::new(1, f64::INFINITY).reciprocal();
        assert_that(z, Complex.ZERO);
    }

    #[test]
    fn testReciprocalZero() {
        assert_that(Complex.ZERO.reciprocal(), Complex.INF);
    }

    #[test]
    fn testReciprocalNaN() {
        assert!(Complex.NaN.reciprocal().isNaN());
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
    fn testMultiplyNaN() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.multiply(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = Complex.NaN.multiply(5);
        Assert.assertSame(Complex.NaN, z);
    }

    #[test]
    fn testMultiplyInfInf() {
        // assert!(infInf.multiply(infInf).isNaN()); // MATH-620
        assert!(infInf.multiply(infInf).isInfinite());
    }

    #[test]
    fn testMultiplyNaNInf() {
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

        w = negInfNegInf.multiply(oneNaN);
        assert!(f64::is_nan(w.get_real()));
        assert!(f64::is_nan(w.get_imaginary()));

        z = Complex::new(1, f64::INFINITY);
        Assert.assertSame(Complex.INF, z.multiply(z));
    }

    #[test]
    fn testScalarMultiply() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = 2.0;
        let yComplex = Complex::new(yDouble);
        assert_that(x.multiply(yComplex), x.multiply(yDouble));
        zInt = -5;
        let zComplex = Complex::new(zInt);
        assert_that(x.multiply(zComplex), x.multiply(zInt));
    }

    #[test]
    fn testScalarMultiplyNaN() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = Double.NaN;
        let yComplex = Complex::new(yDouble);
        assert_that(x.multiply(yComplex), x.multiply(yDouble));
    }

    #[test]
    fn testScalarMultiplyInf() {
        let x = Complex::new(1, 1);
        let yDouble: f64 = Double.POSITIVE_INFINITY;
        let yComplex = Complex::new(yDouble);
        assert_that(x.multiply(yComplex), x.multiply(yDouble));

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = Complex::new(yDouble);
        assert_that(x.multiply(yComplex), x.multiply(yDouble));
    }

    #[test]
    fn testNegate() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.negate();
        assert_that(-3.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(-4.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testNegateNaN() {
        let z: Complex = Complex.NaN.negate();
        assert!(z.isNaN());
    }

    #[test]
    fn testSubtract() {
        let x = Complex::new(3.0, 4.0);
        let y = Complex::new(5.0, 6.0);
        let z: Complex = x.subtract(y);
        assert_that(-2.0, is(close_to(z.get_real(), 1.0e-5)));
        assert_that(-2.0, is(close_to(z.get_imaginary(), 1.0e-5)));
    }

    #[test]
    fn testSubtractNaN() {
        let x = Complex::new(3.0, 4.0);
        let z: Complex = x.subtract(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = Complex::new(1, f64::NAN);
        let w: Complex = x.subtract(z);
        Assert.assertSame(Complex.NaN, w);
    }

    #[test]
    fn testSubtractInf() {
        let x = Complex::new(1, 1);
        let z = Complex::new(f64::INFINITY, 0);
        let w: Complex = x.subtract(z);
        assert_that(w.get_imaginary(), is(close_to(1, 0)));
        assert_that(f64::INFINITY, is(close_to(w.get_real(), 0)));

        x = Complex::new(f64::INFINITY, 0);
        assert!(f64::is_nan(x.subtract(z).get_real()));
    }

    #[test]
    fn testScalarSubtract() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = 2.0;
        let yComplex = Complex::new(yDouble);
        assert_that(x.subtract(yComplex), x.subtract(yDouble));
    }

    #[test]
    fn testScalarSubtractNaN() {
        let x = Complex::new(3.0, 4.0);
        let yDouble: f64 = Double.NaN;
        let yComplex = Complex::new(yDouble);
        assert_that(x.subtract(yComplex), x.subtract(yDouble));
    }

    #[test]
    fn testScalarSubtractInf() {
        let x = Complex::new(1, 1);
        let yDouble: f64 = f64::INFINITY;
        let yComplex = Complex::new(yDouble);
        assert_that(x.subtract(yComplex), x.subtract(yDouble));

        x = Complex::new(f64::INFINITY, 0);
        assert_that(x.subtract(yComplex), x.subtract(yDouble));
    }
}


