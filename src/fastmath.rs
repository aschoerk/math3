#![allow(dead_code)] 
#![allow(unused_parens)]
use std::f64;
use std::f32;
use rsutils::*;
use precision;

/** Archimede's constant PI, ratio of circle circumference to diameter. */
pub const PI: f64 = 105414357.0 / 33554432.0 + 1.984187159361080883e-9;
/** Napier's constant e, base of the natural logarithm. */
pub const E: f64 = 2850325.0 / 1048576.0 + 8.254840070411028747e-8;

/** Index of exp(0) in the array of integer exponentials. */
const EXP_INT_TABLE_MAX_INDEX: i16 = 750;
/** Length of the array of integer exponentials. */
const EXP_INT_TABLE_LEN: usize = (EXP_INT_TABLE_MAX_INDEX * 2) as usize;
/** Logarithm table length. */
const LN_MANT_LEN: usize = 1024;
/** Exponential fractions table length. */
const EXP_FRAC_TABLE_LEN: usize = 1025; // 0, 1/1024, ... 1024/1024

/** StrictMath.log(Double.MAX_VALUE): {@value} */
// const LOG_MAX_VALUE: f64 = f64::MAX.ln();

/** Indicator for tables initialization.
 * <p>
 * This compile-time constant should be set to true only if one explicitly
 * wants to compute the tables at class loading time instead of using the
 * already computed ones provided as literal arrays below.
 * </p>
 */
const RECOMPUTE_TABLES_AT_RUNTIME: bool = false;

/** log(2) (high bits). */
const LN_2_A: f64 = 0.693147063255310059;

/** log(2) (low bits). */
const LN_2_B: f64 = 1.17304635250823482e-7;

 /** Coefficients for log, when input 0.99 < x < 1.01. */
const LN_QUICK_COEF: [[f64; 2]; 9] = [
        [1.0, 5.669184079525E-24],
        [-0.25, -0.25],
        [0.3333333134651184, 1.986821492305628E-8],
        [-0.25, -6.663542893624021E-14],
        [0.19999998807907104, 1.1921056801463227E-8],
        [-0.1666666567325592, -7.800414592973399E-9],
        [0.1428571343421936, 5.650007086920087E-9],
        [-0.12502530217170715, -7.44321345601866E-11],
        [0.11113807559013367, 9.219544613762692E-9],
    ];

    /** Coefficients for log in the range of 1.0 < x < 1.0 + 2^-10. */
const LN_HI_PREC_COEF: [[f64;2];6] = [
    [1.0, -6.032174644509064E-23],
    [-0.25, -0.25],
    [0.3333333134651184, 1.9868161777724352E-8],
    [-0.2499999701976776, -2.957007209750105E-8],
    [0.19999954104423523, 1.5830993332061267E-10],
    [-0.16624879837036133, -2.6033824355191673E-8]
];

/** Sine, Cosine, Tangent tables are for 0, 1/8, 2/8, ... 13/8 = PI/2 approx. */
    const SINE_TABLE_LEN: u16 = 14;

    /** Sine table (high bits). */
    const SINE_TABLE_A: [f64; 14] =
        [
         0.0,
         0.1246747374534607,
         0.24740394949913025,
         0.366272509098053,
         0.4794255495071411,
         0.5850973129272461,
         0.6816387176513672,
         0.7675435543060303,
         0.8414709568023682,
         0.902267575263977,
         0.9489846229553223,
         0.9808930158615112,
         0.9974949359893799,
         0.9985313415527344,
    ];

    /** Sine table (low bits). */
    const SINE_TABLE_B: [f64; 14] =
        [
         0.0,
        -4.068233003401932E-9,
         9.755392680573412E-9,
         1.9987994582857286E-8,
        -1.0902938113007961E-8,
        -3.9986783938944604E-8,
         4.23719669792332E-8,
        -5.207000323380292E-8,
         2.800552834259E-8,
         1.883511811213715E-8,
        -3.5997360512765566E-9,
         4.116164446561962E-8,
         5.0614674548127384E-8,
        -1.0129027912496858E-9,
    ];

    /** Cosine table (high bits). */
    const COSINE_TABLE_A: [f64; 14] =
        [
         1.0,
         0.9921976327896118,
         0.9689123630523682,
         0.9305076599121094,
         0.8775825500488281,
         0.8109631538391113,
         0.7316888570785522,
         0.6409968137741089,
         0.5403022766113281,
         0.4311765432357788,
         0.3153223395347595,
         0.19454771280288696,
         0.07073719799518585,
        -0.05417713522911072,
    ];

    /** Cosine table (low bits). */
    const COSINE_TABLE_B: [f64; 14] =
        [
         0.0,
         3.4439717236742845E-8,
         5.865827662008209E-8,
        -3.7999795083850525E-8,
         1.184154459111628E-8,
        -3.43338934259355E-8,
         1.1795268640216787E-8,
         4.438921624363781E-8,
         2.925681159240093E-8,
        -2.6437112632041807E-8,
         2.2860509143963117E-8,
        -4.813899778443457E-9,
         3.6725170580355583E-9,
         2.0217439756338078E-10,
    ];


    /** Tangent table, used by atan() (high bits). */
    const TANGENT_TABLE_A: [f64; 14] =
        [
         0.0,
         0.1256551444530487,
         0.25534194707870483,
         0.3936265707015991,
         0.5463024377822876,
         0.7214844226837158,
         0.9315965175628662,
         1.1974215507507324,
         1.5574076175689697,
         2.092571258544922,
         3.0095696449279785,
         5.041914939880371,
         14.101419448852539,
        -18.430862426757812,
    ];

    /** Tangent table, used by atan() (low bits). */
    const TANGENT_TABLE_B: [f64; 14] =
        [
         0.0,
        -7.877917738262007E-9,
        -2.5857668567479893E-8,
         5.2240336371356666E-9,
         5.206150291559893E-8,
         1.8307188599677033E-8,
        -5.7618793749770706E-8,
         7.848361555046424E-8,
         1.0708593250394448E-7,
         1.7827257129423813E-8,
         2.893485277253286E-8,
         3.1660099222737955E-7,
         4.983191803254889E-7,
        -3.356118100840571E-7,
    ];
    
       /** Bits of 1/(2*pi), need for reducePayneHanek(). */
    const  RECIP_2PI: [u64; 18] = [
        (0x28be60db << 32) | 0x9391054a,
        (0x7f09d5f4 << 32) | 0x7d4d3770,
        (0x36d8a566 << 32) | 0x4f10e410,
        (0x7f9458ea << 32) | 0xf7aef158,
        (0x6dc91b8e << 32) | 0x909374b8,
        (0x01924bba << 32) | 0x82746487,
        (0x3f877ac7 << 32) | 0x2c4a69cf,
        (0xba208d7d << 32) | 0x4baed121,
        (0x3a671c09 << 32) | 0xad17df90,
        (0x4e64758e << 32) | 0x60d4ce7d,
        (0x272117e2 << 32) | 0xef7e4a0e,
        (0xc7fe25ff << 32) | 0xf7816603,
        (0xfbcbc462 << 32) | 0xd6829b47,
        (0xdb4d9fb3 << 32) | 0xc9f2c26d,
        (0xd3d18fd9 << 32) | 0xa797fa8b,
        (0x5d49eeb1 << 32) | 0xfaf97c5e,
        (0xcf41ce7d << 32) | 0xe294a4ba,
         0x9afed7ec << 32  ];

    /** Bits of pi/4, need for reducePayneHanek(). */
    const PI_O_4_BITS: [u64; 2] = [
        (0xc90fdaa2 << 32) | 0x2168c234,
        (0xc4c6628b << 32) | 0x80dc1cd1 ];

    /** Eighths.
     * This is used by sinQ, because its faster to do a table lookup than
     * a multiply in this time-critical routine
     */
    const EIGHTHS: [f64; 14] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625];

    /** Table of 2^((n+2)/3) */
    const CBRTTWO: [f64; 5] = [ 0.6299605249474366,
                                            0.7937005259840998,
                                            1.0,
                                            1.2599210498948732,
                                            1.5874010519681994 ];

    /*
     *  There are 52 bits in the mantissa of a double.
     *  For additional precision, the code splits double numbers into two parts,
     *  by clearing the low order 30 bits if possible, and then performs the arithmetic
     *  on each half separately.
     */

    /**
     * 0x40000000 - used to split a double into two parts, both with the low order bits cleared.
     * Equivalent to 2^30.
     */
    const HEX_40000000: u64 = 0x40000000; // 1073741824L

    /** Mask used to clear low order 30 bits */
    const MASK_30BITS: u64 = 0xFFFFFFFFFFFFFFFF - (HEX_40000000 -1); // 0xFFFFFFFFC0000000L;

    /** Mask used to clear the non-sign part of an int. */
    const MASK_NON_SIGN_INT: u32 = 0x7fffffff;

    /** Mask used to clear the non-sign part of a long. */
    const MASK_NON_SIGN_LONG: u64 = 0x7fffffffffffffff;

    /** Mask used to extract exponent from double bits. */
    const MASK_DOUBLE_EXPONENT: u64 = 0x7ff0000000000000;

    /** Mask used to extract mantissa from double bits. */
    const MASK_DOUBLE_MANTISSA: u64 = 0x000fffffffffffff;

    /** Mask used to add implicit high order bit for normalized double. */
    const IMPLICIT_HIGH_BIT: u64 = 0x0010000000000000;

    /** 2^52 - double numbers this large must be integral (no fraction) or NaN or Infinite */
    const TWO_POWER_52: f64 = 4503599627370496.0;

    /** Constant: {@value}. */
    const F_1_3: f64 = 1.0 / 3.0;
    /** Constant: {@value}. */
    const F_1_5: f64 = 1.0 / 5.0;
    /** Constant: {@value}. */
    const F_1_7: f64 = 1.0 / 7.0;
    /** Constant: {@value}. */
    const F_1_9: f64 = 1.0 / 9.0;
    /** Constant: {@value}. */
    const F_1_11: f64 = 1.0 / 11.0;
    /** Constant: {@value}. */
    const F_1_13: f64 = 1.0 / 13.0;
    /** Constant: {@value}. */
    const F_1_15: f64 = 1.0 / 15.0;
    /** Constant: {@value}. */
    const F_1_17: f64 = 1.0 / 17.0;
    /** Constant: {@value}. */
    const F_3_4: f64 = 3.0 / 4.0;
    /** Constant: {@value}. */
    const F_15_16: f64 = 15.0 / 16.0;
    /** Constant: {@value}. */
    const F_13_14: f64 = 13.0 / 14.0;
    /** Constant: {@value}. */
    const F_11_12: f64 = 11.0 / 12.0;
    /** Constant: {@value}. */
    const F_9_10: f64 = 9.0 / 10.0;
    /** Constant: {@value}. */
    const F_7_8: f64 = 7.0 / 8.0;
    /** Constant: {@value}. */
    const F_5_6: f64 = 5.0 / 6.0;
    /** Constant: {@value}. */
    const F_1_2: f64 = 1.0 / 2.0;
    /** Constant: {@value}. */
    const F_1_4: f64 = 1.0 / 4.0;
    
    fn a(i: i32) -> [f64; EXP_INT_TABLE_LEN] {
    	[0.0;EXP_INT_TABLE_LEN]
    }
    fn b(i: i32) -> [f64; EXP_FRAC_TABLE_LEN] {
    	[0.0;EXP_FRAC_TABLE_LEN]
    }
    
    lazy_static! {
    	pub static ref EXP_INT_TABLE_A: [f64;EXP_INT_TABLE_LEN] = a(1);
    	pub static ref EXP_INT_TABLE_B: [f64;EXP_INT_TABLE_LEN] = a(2);
    	pub static ref EXP_FRAC_TABLE_A: [f64;EXP_FRAC_TABLE_LEN] = b(1);
    	pub static ref EXP_FRAC_TABLE_B: [f64;EXP_FRAC_TABLE_LEN] = b(2);
    	pub static ref LOG_MAX_VALUE: f64 = f64::ln(f64::MAX);
    }

pub struct F64 {
	
} 
impl F64 {
	 /**
     * Get the high order bits from the mantissa.
     * Equivalent to adding and subtracting HEX_40000 but also works for very large numbers
     *
     * @param d the value to split
     * @return the high order part of the mantissa
     */
    fn double_high_part(d: f64) -> f64 {
        if (d > -*precision::SAFE_MIN && d < *precision::SAFE_MIN){
            return d; // These are un-normalised - don't try to convert
        }
        let xl = double_to_raw_long_bits(&d); // can take raw bits because just gonna convert it back
        xl &= MASK_30BITS; // Drop low order bits
        long_bits_to_double(&xl)
    }
    
     /** Compute the square root of a number.
     * <p><b>Note:</b> this implementation currently delegates to {@link Math#sqrt}
     * @param a number on which evaluation is done
     * @return square root of a
     */
    pub fn sqrt(a: f64) -> f64 {
        f64::sqrt(a)
    }
    
    /**
     * Internal helper method for exponential function.
     * @param x original argument of the exponential function
     * @param extra extra bits of precision on input (To Be Confirmed)
     * @param hiPrec extra bits of precision on output (To Be Confirmed)
     * @return exp(x)
     */
    fn exp_prec(x: f64, extra: f64, hi_prec_o: &mut Option<&mut [f64;2]>) -> f64 {
        let mut int_part_a: f64 = 0.0;
        let mut int_part_b: f64 = 0.0;
        let int_val = x as i16;

        /* Lookup exp(floor(x)).
         * int_part_a will have the upper 22 bits, int_part_b will have the lower
         * 52 bits.
         */
        if x < 0.0 {

            // We don't check against int_val here as conversion of large negative double values
            // may be affected by a JIT bug. Subsequent comparisons can safely use int_val
            if x < -746.0 {
            	match *hi_prec_o {
            		Some (ref mut hi_prec) => {
            			hi_prec[0] = 0.0;
            			hi_prec[1] = 0.0;
            		} ,
            		None => {},
            	}
                return 0.0;
            }

            if (int_val < -709) {
                /* This will produce a subnormal output */
                let result = F64::exp_prec(x+40.19140625, extra, hi_prec_o) / 285040095144011776.0;
            	match *hi_prec_o {
            		Some (ref mut hi_prec) => {
            			hi_prec[0] /= 285040095144011776.0;
            			hi_prec[1] /= 285040095144011776.0;
            		} ,
            		None => {},
            	}                
                return result;
            }

            if (int_val == -709) {
                /* exp(1.494140625) is nearly a machine number... */
                let result = F64::exp_prec(x+1.494140625, extra, hi_prec_o) / 4.455505956692756620;
                match *hi_prec_o {
            		Some (ref mut hi_prec) => {
            			hi_prec[0] /= 4.455505956692756620;
            			hi_prec[1] /= 4.455505956692756620;
            		} ,
            		None => {},
            	}
                return result;
            }

            int_val -= 1;

        } else {
            if (int_val > 709) {
            	match *hi_prec_o {
            		Some (ref mut hi_prec) => {
            			hi_prec[0] = f64::INFINITY;
            			hi_prec[1] = 0.0;
            		} ,
            		None => {},
            	}                
                return f64::INFINITY;
            }

        }

        int_part_a = EXP_INT_TABLE_A[(EXP_INT_TABLE_MAX_INDEX + int_val) as usize];
        int_part_b = EXP_INT_TABLE_B[(EXP_INT_TABLE_MAX_INDEX + int_val) as usize];

        /* Get the fractional part of x, find the greatest multiple of 2^-10 less than
         * x and look up the exp function of it.
         * frac_part_a will have the upper 22 bits, frac_part_b the lower 52 bits.
         */
        let int_frac = ((x - int_val as f64) * 1024.0) as usize;
        let frac_part_a: f64 = EXP_FRAC_TABLE_A[int_frac];
        let frac_part_b: f64 = EXP_FRAC_TABLE_B[int_frac];

        /* epsilon is the difference in x from the nearest multiple of 2^-10.  It
         * has a value in the range 0 <= epsilon < 2^-10.
         * Do the subtraction from x as the last step to avoid possible loss of precision.
         */
        let epsilon: f64 = x - (int_val + int_frac as f64 / 1024.0);

        /* Compute z = exp(epsilon) - 1.0 via a minimax polynomial.  z has
       full double precision (52 bits).  Since z < 2^-10, we will have
       62 bits of precision when combined with the constant 1.  This will be
       used in the last addition below to get proper rounding. */

        /* Remez generated polynomial.  Converges on the interval [0, 2^-10], error
       is less than 0.5 ULP */
        let z: f64 = 0.04168701738764507;
        z = z * epsilon + 0.1666666505023083;
        z = z * epsilon + 0.5000000000042687;
        z = z * epsilon + 1.0;
        z = z * epsilon + -3.940510424527919E-20;

        /* Compute (int_part_a+int_part_b) * (frac_part_a+frac_part_b) by binomial
       expansion.
       tempA is exact since int_part_a and int_part_b only have 22 bits each.
       tempB will have 52 bits of precision.
         */
        let tempA: f64 = int_part_a * frac_part_a;
        let tempB: f64 = int_part_a * frac_part_b + int_part_b * frac_part_a + int_part_b * frac_part_b;

        /* Compute the result.  (1+z)(tempA+tempB).  Order of operations is
       important.  For accuracy add by increasing size.  tempA is exact and
       much larger than the others.  If there are extra bits specified from the
       pow() function, use them. */
        let tempC: f64 = tempB + tempA;

        // If tempC is positive infinite, the evaluation below could result in NaN,
        // because z could be negative at the same time.
        if tempC == f64::INFINITY {
            return f64::INFINITY;
        }

        let result: f64 = 
	        if extra != 0.0 {
	            tempC*extra*z + tempC*extra + tempC*z + tempB + tempA
	        } else {
	            tempC*z + tempB + tempA
	        };

		match *hi_prec_o {
    		Some (ref mut hi_prec) => {
    			hi_prec[0] = tempA;
    			hi_prec[1] = tempC*extra*z + tempC*extra + tempC*z + tempB;
    		} ,
    		None => {},
    	}

        result
    }
    
    /**
     * Exponential function.
     *
     * Computes exp(x), function result is nearly rounded.   It will be correctly
     * rounded to the theoretical value for 99.9% of input values, otherwise it will
     * have a 1 ULP error.
     *
     * Method:
     *    Lookup intVal = exp(int(x))
     *    Lookup fracVal = exp(int(x-int(x) / 1024.0) * 1024.0 );
     *    Compute z as the exponential of the remaining bits by a polynomial minus one
     *    exp(x) = intVal * fracVal * (1 + z)
     *
     * Accuracy:
     *    Calculation is done with 63 bits of precision, so result should be correctly
     *    rounded for 99.9% of input values, with less than 1 ULP error otherwise.
     *
     * @param x   a double
     * @return double e<sup>x</sup>
     */
    pub fn exp(x: f64) -> f64 {
        return F64::exp_prec(x, 0.0, None);
    }

    
     /** Compute the hyperbolic cosine of a number.
     * @param x number on which evaluation is done
     * @return hyperbolic cosine of x
     */
    pub fn cosh(x: f64) -> f64 {
      if x != x {
          return x;
      }

      // cosh[z] = (exp(z) + exp(-z))/2

      // for numbers with magnitude 20 or so,
      // exp(-z) can be ignored in comparison with exp(z)

      if x > 20 {
          if x >= LOG_MAX_VALUE {
              // Avoid overflow (MATH-905).
              let t = F64::exp(0.5 * x);
              (0.5 * t) * t
          } else {
               0.5 * F64::exp(x)
          }
      } else if x < -20 {
          if x <= -LOG_MAX_VALUE {
              // Avoid overflow (MATH-905).
              let t = F64::exp(-0.5 * x);
              (0.5 * t) * t
          } else {
              0.5 * F64::exp(-x)
          }
      }

      let hi_prec: [f64; 2] = [0.0, 0.0];
      if (x < 0.0) {
          x = -x;
      }
      F64::exp_prec(x, 0.0, Some(hi_prec));

      let mut ya = hi_prec[0] + hi_prec[1];
      let mut yb = -(ya - hi_prec[0] - hi_prec[1]);

      let mut temp = ya * HEX_40000000;
      let yaa = ya + temp - temp;
      let yab = ya - yaa;

      // recip = 1/y
      let recip = 1.0/ya;
      let mut temp = recip * HEX_40000000;
      let recipa = recip + temp - temp;
      let mut recipb = recip - recipa;

      // Correct for rounding in division
      recipb += (1.0 - yaa*recipa - yaa*recipb - yab*recipa - yab*recipb) * recip;
      // Account for yb
      recipb += -yb * recip * recip;

      // y = y + 1/y
      temp = ya + recipa;
      yb += -(temp - ya - recipa);
      ya = temp;
      temp = ya + recipb;
      yb += -(temp - ya - recipb);
      ya = temp;

      let result = ya + yb;
      result *= 0.5;
      result
    }

	
	pub fn abs(x: f64) -> f64 {
		long_bits_to_double(&(MASK_NON_SIGN_LONG & double_to_raw_long_bits(&x)))
	}
	/**
     * Compute least significant bit (Unit in Last Position) for a number.
     * @param x number from which ulp is requested
     * @return ulp(x)
     */
	pub fn ulp(x: f64) -> f64 {
		if x.is_infinite() {
			f64::INFINITY
		} else {
			F64::abs(x - long_bits_to_double(&(double_to_raw_long_bits(&x) ^ 1)))
		}
	}
	
	 /**
     * Multiply a double number by a power of 2.
     * @param d number to multiply
     * @param n power of 2
     * @return d &times; 2<sup>n</sup>
     */
    pub fn scalb(d: f64, n: i32) -> f64 {

        // first simple and fast handling when 2^n can be represented using normal numbers
        if (n > -1023) && (n < 1024) {
            return d * long_bits_to_double(&(((n + 0x3ff) as u64) << 52));
        } 

        // handle special cases
        if d.is_nan() || d.is_infinite() || (d == 0.0) {
            return d;
        }
        
        if n < -2098 {
            return if d > 0.0 { 0.0 } else { -0.0 };
        }
        if n > 2097 {
            return if d > 0.0 { f64::INFINITY } else {  f64::NEG_INFINITY };
        }

        // decompose d
        let bits = double_to_raw_long_bits(&d);
        let sign = bits & 0x8000000000000000u64;
        let exponent = ((bits as i64 >> 52) as i32) & 0x7ff;
        let mut mantissa   = bits & 0x000fffffffffffff;

        // compute scaled exponent
        let mut scaled_exponent = exponent + n;

        if n < 0 {
            // we are really in the case n <= -1023
            if scaled_exponent > 0 {
                // both the input and the result are normal numbers, we only adjust the exponent
                long_bits_to_double(&(sign | ((scaled_exponent as u64) << 52) | mantissa))
            } else if scaled_exponent > -53 {
                // the input is a normal number and the result is a subnormal number

                // recover the hidden mantissa bit
                mantissa = mantissa | (1u64 << 52);

                // scales down complete mantissa, hence losing least significant bits
                let most_significant_lost_bit = mantissa & (1u64 << (-scaled_exponent));
                mantissa >>= 1 - scaled_exponent;
                if most_significant_lost_bit != 0 {
                    // we need to add 1 bit to round up the result
                    mantissa += 1;
                }
                long_bits_to_double(&(sign | mantissa))

            } else {
                // no need to compute the mantissa, the number scales down to 0
                if sign == 0 { 0.0 } else { -0.0 }
            }
        } else {
            // we are really in the case n >= 1024
            if exponent == 0 {

                // the input number is subnormal, normalize it
                while (mantissa >> 52) != 1 {
                    mantissa <<= 1;
                    --scaled_exponent;
                }
                scaled_exponent += 1;
                
                mantissa = mantissa & 0x000fffffffffffff;

                if scaled_exponent < 2047 {
                    long_bits_to_double(&(sign | ((scaled_exponent as u64) << 52) | mantissa))
                } else {
                    if sign == 0 { f64::INFINITY } else { f64::NEG_INFINITY }
                }

            } else if scaled_exponent < 2047 {
                long_bits_to_double(&(sign | ((scaled_exponent as u64) << 52) | mantissa))
            } else {
                if sign == 0 { f64::INFINITY } else { f64::NEG_INFINITY }
            }
        }
    }  // scalb
    
    
    /**
     * Get the next machine representable number after a number, moving
     * in the direction of another number.
     * <p>
     * The ordering is as follows (increasing):
     * <ul>
     * <li>-INFINITY</li>
     * <li>-MAX_VALUE</li>
     * <li>-MIN_VALUE</li>
     * <li>-0.0</li>
     * <li>+0.0</li>
     * <li>+MIN_VALUE</li>
     * <li>+MAX_VALUE</li>
     * <li>+INFINITY</li>
     * <li></li>
     * <p>
     * If arguments compare equal, then the second argument is returned.
     * <p>
     * If {@code direction} is greater than {@code d},
     * the smallest machine representable number strictly greater than
     * {@code d} is returned; if less, then the largest representable number
     * strictly less than {@code d} is returned.</p>
     * <p>
     * If {@code d} is infinite and direction does not
     * bring it back to finite numbers, it is returned unchanged.</p>
     *
     * @param d base number
     * @param direction (the only important thing is whether
     * {@code direction} is greater or smaller than {@code d})
     * @return the next machine representable number in the specified direction
     */
    pub fn nextAfter(d: f64, direction: f64) -> f64 {
        // handling of some important special cases
        if (d.is_nan() || direction.is_nan()) {
            return f64::NAN;
        } else if (d == direction) {
            return direction;
        } else if (d.is_infinite()) {
            return if d < 0.0 { -f64::MAX } else  { f64::MAX };
        } else if (d == 0.0) {
            return if direction < 0.0 { -f64::MIN_POSITIVE } else { f64::MIN_POSITIVE };
        }
        // special cases MAX_VALUE to infinity and  MIN_VALUE to 0
        // are handled just as normal numbers
        // can use raw bits since already dealt with infinity and NaN
        let bits = double_to_raw_long_bits(&d);
        let sign = bits & 0x8000000000000000u64;
        if ((direction < d) ^ (sign == 0)) {
            return long_bits_to_double(sign | ((bits & 0x7fffffffffffffff) + 1));
        } else {
            return long_bits_to_double(sign | ((bits & 0x7fffffffffffffff) - 1));
        }   
    }
    
    
}

pub struct F32 {
	
	
} 
impl F32 {
	pub fn abs(x: f32) -> f32 {
		int_bits_to_float(&(MASK_NON_SIGN_INT & float_to_raw_int_bits(&x)))
	}
	/**
     * Compute least significant bit (Unit in Last Position) for a number.
     * @param x number from which ulp is requested
     * @return ulp(x)
     */
	pub fn ulp(x: f32) -> f32 {
		if x.is_infinite() {
			f32::INFINITY
		} else {
			F32::abs(x - int_bits_to_float(&(float_to_raw_int_bits(&x) ^ 1)))
		}
	}
	
	 /**
     * Multiply a double number by a power of 2.
     * @param d number to multiply
     * @param n power of 2
     * @return d &times; 2<sup>n</sup>
     */
    pub fn scalb(d: f32, n: i32) -> f32 {

        // first simple and fast handling when 2^n can be represented using normal numbers
        if (n > -127) && (n < 128) {
            return d * int_bits_to_float(&(((n + 127) as u32) << 23));
        } 

        // handle special cases
        if d.is_nan() || d.is_infinite() || (d == 0.0) {
            return d;
        }
        
        if n < -277 {
            return if d > 0.0 { 0.0 } else { -0.0 };
        }
        if n > 276 {
            return if d > 0.0 { f32::INFINITY } else {  f32::NEG_INFINITY };
        }

        // decompose d
        let bits = double_to_raw_long_bits(&d);
        let sign = bits & 0x80000000u32;
        let exponent = ((bits as i32 >> 23) as i32) & 0xff;
        let mut mantissa   = bits & 0x007fffff;

        // compute scaled exponent
        let mut scaled_exponent = exponent + n;

        if n < 0 {
            // we are really in the case n <= -127
            if scaled_exponent > 0 {
                // both the input and the result are normal numbers, we only adjust the exponent
                long_bits_to_double(&(sign | ((scaled_exponent as u32) << 23) | mantissa))
            } else if scaled_exponent > -24 {
                // the input is a normal number and the result is a subnormal number

                // recover the hidden mantissa bit
                mantissa = mantissa | (1u32 << 23);

                // scales down complete mantissa, hence losing least significant bits
                let most_significant_lost_bit = mantissa & (1u32 << (-scaled_exponent));
                mantissa >>= 1 - scaled_exponent;
                if most_significant_lost_bit != 0 {
                    // we need to add 1 bit to round up the result
                    mantissa += 1;
                }
                long_bits_to_double(&(sign | mantissa))

            } else {
                // no need to compute the mantissa, the number scales down to 0
                if sign == 0 { 0.0 } else { -0.0 }
            }
        } else {
            // we are really in the case n >= 128
            if exponent == 0 {

                // the input number is subnormal, normalize it
                while (mantissa >> 23) != 1 {
                    mantissa <<= 1;
                    --scaled_exponent;
                }
                scaled_exponent += 1;
                
                mantissa = mantissa & 0x007fffff;

                if scaled_exponent < 255 {
                    long_bits_to_double(&(sign | ((scaled_exponent as u32) << 23) | mantissa))
                } else {
                    if sign == 0 { f32::INFINITY } else { f32::NEG_INFINITY }
                }

            } else if scaled_exponent < 255 {
                long_bits_to_double(&(sign | ((scaled_exponent as u32) << 23) | mantissa))
            } else {
                if sign == 0 { f32::INFINITY } else { f32::NEG_INFINITY }
            }
        }
    }  // scalb
    
        /**
     * Get the next machine representable number after a number, moving
     * in the direction of another number.
     * <p>
     * The ordering is as follows (increasing):
     * <ul>
     * <li>-INFINITY</li>
     * <li>-MAX_VALUE</li>
     * <li>-MIN_VALUE</li>
     * <li>-0.0</li>
     * <li>+0.0</li>
     * <li>+MIN_VALUE</li>
     * <li>+MAX_VALUE</li>
     * <li>+INFINITY</li>
     * <li></li>
     * <p>
     * If arguments compare equal, then the second argument is returned.
     * <p>
     * If {@code direction} is greater than {@code f},
     * the smallest machine representable number strictly greater than
     * {@code f} is returned; if less, then the largest representable number
     * strictly less than {@code f} is returned.</p>
     * <p>
     * If {@code f} is infinite and direction does not
     * bring it back to finite numbers, it is returned unchanged.</p>
     *
     * @param f base number
     * @param direction (the only important thing is whether
     * {@code direction} is greater or smaller than {@code f})
     * @return the next machine representable number in the specified direction
     */
     pub fn nextAfter(d: f32, direction: f32) -> f32 {
		// handling of some important special cases
        if (d.is_nan() || direction.is_nan()) {
            return f32::NAN;
        } else if (d == direction) {
            return direction;
        } else if (d.is_infinite()) {
            return if d < 0.0 { -f32::MAX } else  { f32::MAX };
        } else if (d == 0.0) {
            return if direction < 0.0 { -f32::MIN_POSITIVE } else { f32::MIN_POSITIVE };
        }
        // special cases MAX_VALUE to infinity and  MIN_VALUE to 0
        // are handled just as normal numbers
        // can use raw bits since already dealt with infinity and NaN
        let bits = float_to_raw_int_bits(&d);
        let sign = bits & 0x80000000u64;
        if ((direction < d) ^ (sign == 0)) {
            return int_bits_to_float(sign | ((bits & 0x7fffffff) + 1));
        } else {
            return int_bits_to_float(sign | ((bits & 0x7fffffff) - 1));
        }
        
    }
    
}


pub struct I32 {
}
impl I32 {
	pub fn abs(x: i32) -> i32 {
		let i = x >> 31;
		(x ^ ((i as u32 ^ 0xFFFFFFFF) as i32 + 1)) + i
	} 
}

pub struct I64 {
}

impl I64 {
	pub fn abs(x: i64) -> i64 {
		let l = x >> 63;
		(x ^ ((l as u64 ^ 0xFFFFFFFFFFFFFFFF) as i64 + 1)) + l
	} 
}

#[cfg(test)]
mod tests {
	use fastmath::F64;
	use fastmath::F32;
	use fastmath::I64;
	use fastmath::I32;
	
	#[test]
	fn f64_abs() {
		assert!(F64::abs(-1.0) == 1.0);
		assert!(F64::abs(-1.0e7) == 1.0e7);
		assert!(F64::abs(1.0) == 1.0);
		assert!(F64::abs(1.0e7) == 1.0e7);
		assert!(F64::abs(-1.0e-7) == 1.0e-7);
		assert!(F64::abs(-0.0) == 0.0);
		assert!(F64::abs(-1.011) == 1.011);		
	}
	
	#[test]
	fn f32_abs() {
		assert!(F32::abs(-1.0) == 1.0);
		assert!(F32::abs(-1.0e7) == 1.0e7);
		assert!(F32::abs(1.0) == 1.0);
		assert!(F32::abs(1.0e7) == 1.0e7);
		assert!(F32::abs(-1.0e-7) == 1.0e-7);
		assert!(F32::abs(-0.0) == 0.0);
		assert!(F32::abs(-1.011) == 1.011);		
	}
}
