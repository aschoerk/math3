use base::*;    
use dfp::dfp::*;            
                   
/** IEEE 854-1987 flag for invalid operation. */
 const FLAG_INVALID: i32 = 1;

/** IEEE 854-1987 flag for division by zero. */
 const FLAG_DIV_ZERO: i32 = 2;

/** IEEE 854-1987 flag for overflow. */
 const FLAG_OVERFLOW: i32 = 4;

/** IEEE 854-1987 flag for underflow. */
 const FLAG_UNDERFLOW: i32 = 8;

/** IEEE 854-1987 flag for inexact result. */
 const FLAG_INEXACT: i32 = 16;

lazy_static! {
	
//    pub static ref statics =  {
//        if sqr2_string == null || sqr2_string.length() < high_precision_decimal_digits - 3 {
//        	let high_precision_decimal_digits: i32 = 200;
//            // recompute the string representation of the transcendental constants
//             let high_precision_field: DfpField = Dfp::newField(high_precision_decimal_digits, false);
//             let high_precision_one: Dfp = Dfp::new(high_precision_field, 1);
//             let high_precision_two: Dfp = Dfp::new(high_precision_field, 2);
//             let high_precision_three: Dfp = Dfp::new(high_precision_field, 3);
//             let high_precision_sqr2: Dfp = high_precision_two.sqrt();
//             let high_precision_sqr3: Dfp = high_precision_three.sqrt();
//            (high_precision_sqr2.to_string(),
//            high_precision_one.divide(high_precision_sqr2).to_string(),
//            high_precision_sqr3.to_string(),
//            high_precision_one.divide(high_precision_sqr3).to_string(),
//            compute_pi(high_precision_one, high_precision_two, high_precision_three).to_string(),
//            compute_exp(high_precision_one, high_precision_one).to_string(),
//            compute_ln(high_precision_two, high_precision_one, high_precision_two).to_string(),
//            compute_ln(Dfp::new(high_precision_field, 5), high_precision_one, high_precision_two).to_string(),
//            compute_ln(Dfp::new(high_precision_field, 10), high_precision_one, high_precision_two).to_string())
//        }
//    }
//	/** High precision string representation of √2. */
	 pub static ref sqr2_string: String = "".to_string();
	
	// Note: the static strings are set up (once) by the ctor and @GuardedBy("DfpField.class")
	/** High precision string representation of √2 / 2. */
	 pub static ref sqr2_reciprocal_string: String = "".to_string();
	
	/** High precision string representation of √3. */
	pub static ref sqr3_string: String = "".to_string();
	
	/** High precision string representation of √3 / 3. */
	pub static ref sqr3_reciprocal_string: String = "".to_string();
	
	/** High precision string representation of π. */
	pub static ref pi_string: String = "".to_string();
	
	/** High precision string representation of e. */
	pub static ref e_string: String = "".to_string();
	
	/** High precision string representation of ln(2). */
	pub static ref ln2_string: String = "".to_string();
	
	/** High precision string representation of ln(5). */
	pub static ref ln5_string: String = "".to_string();
	
	/** High precision string representation of ln(10). */
	pub static ref ln10_string: String = "".to_string();
}


pub struct DfpField {

    /** The number of radix digits.
     * Note these depend on the radix which is 10000 digits,
     * so each one is equivalent to 4 decimal digits.
     */
     radix_digits: i32, 

    /** A {@link Dfp} with value 0. */
     zero: Dfp,

    /** A {@link Dfp} with value 1. */
     one: Dfp,

    /** A {@link Dfp} with value 2. */
     two: Dfp,

    /** A {@link Dfp} with value √2. */
     sqr2: Dfp,

    /** A two elements {@link Dfp} array with value √2 split in two pieces. */
     sqr2_split: [Dfp; 2],

    /** A {@link Dfp} with value √2 / 2. */
     sqr2_reciprocal: Dfp,

    /** A {@link Dfp} with value √3. */
     sqr3: Dfp,

    /** A {@link Dfp} with value √3 / 3. */
     sqr3_reciprocal: Dfp,

    /** A {@link Dfp} with value π. */
     pi: Dfp,

    /** A two elements {@link Dfp} array with value π split in two pieces. */
     pi_split: [Dfp; 2],

    /** A {@link Dfp} with value e. */
     e: Dfp,

    /** A two elements {@link Dfp} array with value e split in two pieces. */
     e_split: [Dfp; 2],

    /** A {@link Dfp} with value ln(2). */
     ln2: Dfp,

    /** A two elements {@link Dfp} array with value ln(2) split in two pieces. */
     ln2_split: [Dfp; 2],

    /** A {@link Dfp} with value ln(5). */
     ln5: Dfp,

    /** A two elements {@link Dfp} array with value ln(5) split in two pieces. */
     ln5_split: [Dfp; 2],

    /** A {@link Dfp} with value ln(10). */
     ln10: Dfp,

    /** Current rounding mode. */
     r_mode: RoundingMode,

    /** IEEE 854-1987 signals. */
     ieee_flags: i32,
}

/** Enumerate for rounding modes. */
pub enum RoundingMode {

    /** Rounds toward zero (truncation). */
    ROUND_DOWN, /** Rounds away from zero if discarded digit is non-zero. */
    ROUND_UP, /** Rounds towards nearest unless both are equidistant in which case it rounds away from zero. */
    ROUND_HALF_UP, /** Rounds towards nearest unless both are equidistant in which case it rounds toward zero. */
    ROUND_HALF_DOWN, /** Rounds towards nearest unless both are equidistant in which case it rounds toward the even neighbor.
     * This is the default as  specified by IEEE 854-1987
     */
    ROUND_HALF_EVEN, /** Rounds towards nearest unless both are equidistant in which case it rounds toward the odd neighbor.  */
    ROUND_HALF_ODD, /** Rounds towards positive infinity. */
    ROUND_CEIL, /** Rounds towards negative infinity. */
    ROUND_FLOOR
}


impl DfpField {


    /** Create a factory for the specified number of radix digits.
     * <p>
     * Note that since the {@link Dfp} class uses 10000 as its radix, each radix
     * digit is equivalent to 4 decimal digits. This implies that asking for
     * 13, 14, 15 or 16 decimal digits will really lead to a 4 radix 10000 digits in
     * all cases.
     * </p>
     * @param decimalDigits minimal number of decimal digits.
     */
    pub fn new( decimal_digits: i32) -> DfpField {
        self(decimal_digits, true);
    }

    /** Create a factory for the specified number of radix digits.
     * <p>
     * Note that since the {@link Dfp} class uses 10000 as its radix, each radix
     * digit is equivalent to 4 decimal digits. This implies that asking for
     * 13, 14, 15 or 16 decimal digits will really lead to a 4 radix 10000 digits in
     * all cases.
     * </p>
     * @param decimalDigits minimal number of decimal digits
     * @param computeConstants if true, the transcendental constants for the given precision
     * must be computed (setting this flag to false is RESERVED for the internal recursive call)
     */
    fn new( decimal_digits: i32,  compute_constants: bool) -> DfpField {
        self.radixDigits =  if (decimal_digits < 13) { 4 } else { (decimal_digits + 3) / 4 };
        self.rMode = RoundingMode::ROUND_HALF_EVEN;
        self.ieeeFlags = 0;
        self.zero = Dfp::new_i32(self, 0);
        self.one = Dfp::new_i32(self, 1);
        self.two = Dfp::new_i32(self, 2);
        if compute_constants {
            // set up transcendental constants
         
                // set up the constants at current field accuracy
                sqr2 = Dfp::new_string(self, sqr2_string);
                sqr2_split = split(sqr2_string);
                sqr2_reciprocal = Dfp::new_string(self, sqr2_reciprocal_string);
                sqr3 = Dfp::new(self, sqr3_string);
                sqr3_reciprocal = Dfp::new_string(self, sqr3_reciprocal_string);
                pi = Dfp::new_string(self, pi_string);
                pi_split = split(pi_string);
                e = Dfp::new_string(self, e_string);
                e_split = split(e_string);
                ln2 = Dfp::new_string(self, ln2_string);
                ln2_split = split(ln2_string);
                ln5 = Dfp::new_string(self, ln5_string);
                ln5_split = split(ln5_string);
                ln10 = Dfp::new_string(self, ln10_string);
            
        } else {
            // dummy settings for unused constants
            sqr2 = null;
            sqr2_split = null;
            sqr2_reciprocal = null;
            sqr3 = null;
            sqr3_reciprocal = null;
            pi = null;
            pi_split = null;
            e = null;
            e_split = null;
            ln2 = null;
            ln2_split = null;
            ln5 = null;
            ln5_split = null;
            ln10 = null;
        }
    }

    /** Get the number of radix digits of the {@link Dfp} instances built by self factory.
     * @return number of radix digits
     */
    pub fn  get_radix_digits(&self) -> i32  {
        return self.radix_digits;
    }

    /** Set the rounding mode.
     *  If not set, the default value is {@link RoundingMode#ROUND_HALF_EVEN}.
     * @param mode desired rounding mode
     * Note that the rounding mode is common to all {@link Dfp} instances
     * belonging to the current {@link DfpField} in the system and will
     * affect all future calculations.
     */
    pub fn  set_rounding_mode(&self,  mode: &RoundingMode) -> void  {
        self.r_mode = mode;
    }

    /** Get the current rounding mode.
     * @return current rounding mode
     */
    pub fn  get_rounding_mode(&self) -> RoundingMode  {
        return self.r_mode;
    }

    /** Get the IEEE 854 status flags.
     * @return IEEE 854 status flags
     * @see #clearIEEEFlags()
     * @see #setIEEEFlags(int)
     * @see #setIEEEFlagsBits(int)
     * @see #FLAG_INVALID
     * @see #FLAG_DIV_ZERO
     * @see #FLAG_OVERFLOW
     * @see #FLAG_UNDERFLOW
     * @see #FLAG_INEXACT
     */
    pub fn  get_i_e_e_e_flags(&self) -> i32  {
        return self.ieee_flags;
    }

    /** Clears the IEEE 854 status flags.
     * @see #getIEEEFlags()
     * @see #setIEEEFlags(int)
     * @see #setIEEEFlagsBits(int)
     * @see #FLAG_INVALID
     * @see #FLAG_DIV_ZERO
     * @see #FLAG_OVERFLOW
     * @see #FLAG_UNDERFLOW
     * @see #FLAG_INEXACT
     */
    pub fn  clear_i_e_e_e_flags(&self) -> void  {
        self.ieee_flags = 0;
    }

    /** Sets the IEEE 854 status flags.
     * @param flags desired value for the flags
     * @see #getIEEEFlags()
     * @see #clearIEEEFlags()
     * @see #setIEEEFlagsBits(int)
     * @see #FLAG_INVALID
     * @see #FLAG_DIV_ZERO
     * @see #FLAG_OVERFLOW
     * @see #FLAG_UNDERFLOW
     * @see #FLAG_INEXACT
     */
    pub fn  set_i_e_e_e_flags(&self,  flags: i32) -> void  {
        self.ieee_flags = flags & (FLAG_INVALID | FLAG_DIV_ZERO | FLAG_OVERFLOW | FLAG_UNDERFLOW | FLAG_INEXACT);
    }

    /** Sets some bits in the IEEE 854 status flags, without changing the already set bits.
     * <p>
     * Calling self method is equivalent to call {@code setIEEEFlags(getIEEEFlags() | bits)}
     * </p>
     * @param bits bits to set
     * @see #getIEEEFlags()
     * @see #clearIEEEFlags()
     * @see #setIEEEFlags(int)
     * @see #FLAG_INVALID
     * @see #FLAG_DIV_ZERO
     * @see #FLAG_OVERFLOW
     * @see #FLAG_UNDERFLOW
     * @see #FLAG_INEXACT
     */
    pub fn  set_i_e_e_e_flags_bits(&self,  bits: i32) -> void  {
        self.ieee_flags |= bits & (FLAG_INVALID | FLAG_DIV_ZERO | FLAG_OVERFLOW | FLAG_UNDERFLOW | FLAG_INEXACT);
    }

    /** Makes a {@link Dfp} with a value of 0.
     * @return a new {@link Dfp} with a value of 0
     */
    pub fn  new_dfp(&self) -> Dfp  {
        return Dfp::new(self);
    }

    /** Create an instance from a byte value.
     * @param x value to convert to an instance
     * @return a new {@link Dfp} with the same value as x
     */
    pub fn  new_dfp_i8(&self,  x: i8) -> Dfp  {
        return Dfp::new_i8(self, x);
    }

    /** Create an instance from an int value.
     * @param x value to convert to an instance
     * @return a new {@link Dfp} with the same value as x
     */
    pub fn  new_dfp_i32(&self,  x: i32) -> Dfp  {
        return Dfp::new_i32(self, x);
    }

    /** Create an instance from a long value.
     * @param x value to convert to an instance
     * @return a new {@link Dfp} with the same value as x
     */
    pub fn  new_dfp_i64(&self,  x: i64) -> Dfp  {
        return Dfp::new_i64(self, x);
    }

    /** Create an instance from a double value.
     * @param x value to convert to an instance
     * @return a new {@link Dfp} with the same value as x
     */
    pub fn  new_dfp_f64(&self,  x: f64) -> Dfp  {
        return Dfp::new_f64(self, x);
    }

    /** Copy constructor.
     * @param d instance to copy
     * @return a new {@link Dfp} with the same value as d
     */
    pub fn  new_dfp_dfp(&self,  d: &Dfp) -> Dfp  {
        return Dfp::new_dfp(d);
    }

    /** Create a {@link Dfp} given a String representation.
     * @param s string representation of the instance
     * @return a new {@link Dfp} parsed from specified string
     */
    pub fn  new_dfp_str(&self,  s: &String) -> Dfp  {
        return Dfp::new_str(self, s);
    }

    /** Creates a {@link Dfp} with a non-finite value.
     * @param sign sign of the Dfp to create
     * @param nans code of the value, must be one of {@link Dfp#INFINITE},
     * {@link Dfp#SNAN},  {@link Dfp#QNAN}
     * @return a new {@link Dfp} with a non-finite value
     */
    pub fn  new_dfp(&self,  sign: i8,  nans: i8) -> Dfp  {
        return Dfp::new(self, sign, nans);
    }

    /** Get the constant 0.
     * @return a {@link Dfp} with value 0
     */
    pub fn  get_zero(&self) -> Dfp  {
        return self.zero;
    }

    /** Get the constant 1.
     * @return a {@link Dfp} with value 1
     */
    pub fn  get_one(&self) -> Dfp  {
        return self.one;
    }

    /** {@inheritDoc} */
//    pub fn  get_runtime_class(&self) -> Class<? extends FieldElement<Dfp>>  {
//        return Dfp.class;
//    }

    /** Get the constant 2.
     * @return a {@link Dfp} with value 2
     */
    pub fn  get_two(&self) -> Dfp  {
        return self.two;
    }

    /** Get the constant √2.
     * @return a {@link Dfp} with value √2
     */
    pub fn  get_sqr2(&self) -> Dfp  {
        return self.sqr2;
    }

    /** Get the constant √2 split in two pieces.
     * @return a {@link Dfp} with value √2 split in two pieces
     */
    pub fn  get_sqr2_split(&self) -> [Dfp; 2]  {
        return self.sqr2_split.clone();
    }

    /** Get the constant √2 / 2.
     * @return a {@link Dfp} with value √2 / 2
     */
    pub fn  get_sqr2_reciprocal(&self) -> Dfp  {
        return self.sqr2_reciprocal;
    }

    /** Get the constant √3.
     * @return a {@link Dfp} with value √3
     */
    pub fn  get_sqr3(&self) -> Dfp  {
        return self.sqr3;
    }

    /** Get the constant √3 / 3.
     * @return a {@link Dfp} with value √3 / 3
     */
    pub fn  get_sqr3_reciprocal(&self) -> Dfp  {
        return self.sqr3_reciprocal;
    }

    /** Get the constant π.
     * @return a {@link Dfp} with value π
     */
    pub fn  get_pi(&self) -> Dfp  {
        return self.pi;
    }

    /** Get the constant π split in two pieces.
     * @return a {@link Dfp} with value π split in two pieces
     */
    pub fn  get_pi_split(&self) -> [Dfp; 2]  {
        return self.pi_split.clone();
    }

    /** Get the constant e.
     * @return a {@link Dfp} with value e
     */
    pub fn  get_e(&self) -> Dfp  {
        return self.e;
    }

    /** Get the constant e split in two pieces.
     * @return a {@link Dfp} with value e split in two pieces
     */
    pub fn  get_e_split(&self) -> [Dfp; 2]  {
        return self.e_split.clone();
    }

    /** Get the constant ln(2).
     * @return a {@link Dfp} with value ln(2)
     */
    pub fn  get_ln2(&self) -> Dfp  {
        return self.ln2;
    }

    /** Get the constant ln(2) split in two pieces.
     * @return a {@link Dfp} with value ln(2) split in two pieces
     */
    pub fn  get_ln2_split(&self) -> [Dfp; 2]  {
        return self.ln2_split.clone();
    }

    /** Get the constant ln(5).
     * @return a {@link Dfp} with value ln(5)
     */
    pub fn  get_ln5(&self) -> Dfp  {
        return self.ln5;
    }

    /** Get the constant ln(5) split in two pieces.
     * @return a {@link Dfp} with value ln(5) split in two pieces
     */
    pub fn  get_ln5_split(&self) -> [Dfp; 2]  {
        return self.ln5_split.clone();
    }

    /** Get the constant ln(10).
     * @return a {@link Dfp} with value ln(10)
     */
    pub fn  get_ln10(&self) -> Dfp  {
        return self.ln10;
    }

    /** Breaks a string representation up into two {@link Dfp}'s.
     * The split is such that the sum of them is equivalent to the input string,
     * but has higher precision than using a single Dfp.
     * @param a string representation of the number to split
     * @return an array of two {@link Dfp Dfp} instances which sum equals a
     */
    fn  split(&self,  a: &String) -> [Dfp; 2]  {
         let result: [Dfp; 2];
         let leading: bool = true;
         let sp: i32 = 0;
         let sig: i32 = 0;
         let buf: [char; a.len()];
         for i in 0..a.len() {
            let c = a.char_at(i);
            if c >= '1' && c <= '9' {
                leading = false;
            }
            if c == '.' {
                sig += (400 - sig) % 4;
                leading = false;
            }
            if sig == (radix_digits / 2) * 4 {
                sp = i;
                break;
            }
            if c >= '0' && c <= '9' && !leading {
                sig+= 1;
            }
            buf[i] = c; 
            
        }
        result[0] = Dfp::new_str(self, String::new(buf, 0, sp));
        for i in 0..buf.len() {
        	let c = a.char_at(i);
            buf[i] = Some(c);
            if c >= '0' && c <= '9' && i < sp {
                buf[i] = Some('0');
            }
        }
        result[1] = Dfp::new_str(self, String::new(buf));
        return result;
    }

   

    /** Compute π using Jonathan and Peter Borwein quartic formula.
     * @param one constant with value 1 at desired precision
     * @param two constant with value 2 at desired precision
     * @param three constant with value 3 at desired precision
     * @return π
     */
    fn  compute_pi( one: &Dfp,  two: &Dfp,  three: &Dfp) -> Dfp  {
         let sqrt2: Dfp = two.sqrt();
         let yk: Dfp = sqrt2.subtract(one);
         let four: Dfp = two.add(two);
         let two2kp3: Dfp = two;
         let ak: Dfp = two.multiply(three.subtract(two.multiply(sqrt2)));
        // So the limit here is considered sufficient for most purposes ...
        for i in 1..20 {
             let yk_m1: Dfp = yk;
             let y2: Dfp = yk.multiply(yk);
             let one_minus_y4: Dfp = one.subtract(y2.multiply(y2));
             let s: Dfp = one_minus_y4.sqrt().sqrt();
            yk = one.subtract(s).divide(one.add(s));
            two2kp3 = two2kp3.multiply(four);
             let p: Dfp = one.add(yk);
             let p2: Dfp = p.multiply(p);
            ak = ak.multiply(p2.multiply(p2)).subtract(two2kp3.multiply(yk).multiply(one.add(yk).add(yk.multiply(yk))));
            if yk.equals(yk_m1) {
                break;
            }
        }
        return one.divide(ak);
    }

    /** Compute exp(a).
     * @param a number for which we want the exponential
     * @param one constant with value 1 at desired precision
     * @return exp(a)
     */
    pub fn  compute_exp( a: &Dfp,  one: &Dfp) -> Dfp  {
         let y: Dfp = Dfp::new(one);
         let py: Dfp = Dfp::new(one);
         let f: Dfp = Dfp::new(one);
         let fi: Dfp = Dfp::new(one);
         let x: Dfp = Dfp::new(one);
         for i in 0..10000 {
            x = x.multiply(a);
            y = y.add(x.divide(f));
            fi = fi.add(one);
            f = f.multiply(fi);
            if y.equals(py) {
                break;
            }
            py = Dfp::new(y);
        }
        return y;
    }

    /** Compute ln(a).
     *
     *  Let f(x) = ln(x),
     *
     *  We know that f'(x) = 1/x, thus from Taylor's theorem we have:
     *
     *           -----          n+1         n
     *  f(x) =   \           (-1)    (x - 1)
     *           /          ----------------    for 1 <= n <= infinity
     *           -----             n
     *
     *  or
     *                       2        3       4
     *                   (x-1)   (x-1)    (x-1)
     *  ln(x) =  (x-1) - ----- + ------ - ------ + ...
     *                     2       3        4
     *
     *  alternatively,
     *
     *                  2    3   4
     *                 x    x   x
     *  ln(x+1) =  x - -  + - - - + ...
     *                 2    3   4
     *
     *  This series can be used to compute ln(x), but it converges too slowly.
     *
     *  If we substitute -x for x above, we get
     *
     *                   2    3    4
     *                  x    x    x
     *  ln(1-x) =  -x - -  - -  - - + ...
     *                  2    3    4
     *
     *  Note that all terms are now negative.  Because the even powered ones
     *  absorbed the sign.  Now, subtract the series above from the previous
     *  one to get ln(x+1) - ln(1-x).  Note the even terms cancel out leaving
     *  only the odd ones
     *
     *                             3     5      7
     *                           2x    2x     2x
     *  ln(x+1) - ln(x-1) = 2x + --- + --- + ---- + ...
     *                            3     5      7
     *
     *  By the property of logarithms that ln(a) - ln(b) = ln (a/b) we have:
     *
     *                                3        5        7
     *      x+1           /          x        x        x          \
     *  ln ----- =   2 *  |  x  +   ----  +  ----  +  ---- + ...  |
     *      x-1           \          3        5        7          /
     *
     *  But now we want to find ln(a), so we need to find the value of x
     *  such that a = (x+1)/(x-1).   This is easily solved to find that
     *  x = (a-1)/(a+1).
     * @param a number for which we want the exponential
     * @param one constant with value 1 at desired precision
     * @param two constant with value 2 at desired precision
     * @return ln(a)
     */
    pub fn  compute_ln( a: &Dfp,  one: &Dfp,  two: &Dfp) -> Dfp  {
         let den: i32 = 1;
         let x: Dfp = a.add(Dfp::new(a.get_field(), -1)).divide(a.add(one));
         let y: Dfp = Dfp::new(x);
         let num: Dfp = Dfp::new(x);
         let py: Dfp = Dfp::new(y);
         for i in 0..10000 {
            num = num.multiply(x);
            num = num.multiply(x);
            den += 2;
             let t: Dfp = num.divide(den);
            y = y.add(t);
            if y.equals(py) {
                break;
            }
            py = Dfp::new(y);
        }
        return y.multiply(two);
    }
}


                
                