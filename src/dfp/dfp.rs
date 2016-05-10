                
                   /*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

use base::RealFieldElement;
use util::fastmath;

/**
 *  Decimal floating point library for Java
 *
 *  <p>Another floating point class.  This one is built using radix 10000
 *  which is 10<sup>4</sup>, so its almost decimal.</p>
 *
 *  <p>The design goals here are:
 *  <ol>
 *    <li>Decimal math, or close to it</li>
 *    <li>Settable precision (but no mix between numbers using different settings)</li>
 *    <li>Portability.  Code should be kept as portable as possible.</li>
 *    <li>Performance</li>
 *    <li>Accuracy  - Results should always be +/- 1 ULP for basic
 *         algebraic operation</li>
 *    <li>Comply with IEEE 854-1987 as much as possible.
 *         (See IEEE 854-1987 notes below)</li>
 *  </ol></p>
 *
 *  <p>Trade offs:
 *  <ol>
 *    <li>Memory foot print.  I'm using more memory than necessary to
 *         represent numbers to get better performance.</li>
 *    <li>Digits are bigger, so rounding is a greater loss.  So, if you
 *         really need 12 decimal digits, better use 4 base 10000 digits
 *         there can be one partially filled.</li>
 *  </ol></p>
 *
 *  <p>Numbers are represented  in the following form:
 *  <pre>
 *  n  =  sign × mant × (radix)<sup>exp</sup>;</p>
 *  </pre>
 *  where sign is ±1, mantissa represents a fractional number between
 *  zero and one.  mant[0] is the least significant digit.
 *  exp is in the range of -32767 to 32768</p>
 *
 *  <p>IEEE 854-1987  Notes and differences</p>
 *
 *  <p>IEEE 854 requires the radix to be either 2 or 10.  The radix here is
 *  10000, so that requirement is not met, but  it is possible that a
 *  subclassed can be made to make it behave as a radix 10
 *  number.  It is my opinion that if it looks and behaves as a radix
 *  10 number then it is one and that requirement would be met.</p>
 *
 *  <p>The radix of 10000 was chosen because it should be faster to operate
 *  on 4 decimal digits at once instead of one at a time.  Radix 10 behavior
 *  can be realized by adding an additional rounding step to ensure that
 *  the number of decimal digits represented is constant.</p>
 *
 *  <p>The IEEE standard specifically leaves out internal data encoding,
 *  so it is reasonable to conclude that such a subclass of this radix
 *  10000 system is merely an encoding of a radix 10 system.</p>
 *
 *  <p>IEEE 854 also specifies the existence of "sub-normal" numbers.  This
 *  class does not contain any such entities.  The most significant radix
 *  10000 digit is always non-zero.  Instead, we support "gradual underflow"
 *  by raising the underflow flag for numbers less with exponent less than
 *  expMin, but don't flush to zero until the exponent reaches MIN_EXP-digits.
 *  Thus the smallest number we can represent would be:
 *  1E(-(MIN_EXP-digits-1)*4),  eg, for digits=5, MIN_EXP=-32767, that would
 *  be 1e-131092.</p>
 *
 *  <p>IEEE 854 defines that the implied radix point lies just to the right
 *  of the most significant digit and to the left of the remaining digits.
 *  This implementation puts the implied radix point to the left of all
 *  digits including the most significant one.  The most significant digit
 *  here is the one just to the right of the radix point.  This is a fine
 *  detail and is really only a matter of definition.  Any side effects of
 *  this can be rendered invisible by a subclass.</p>
 * @see DfpField
 * @since 2.2
 */

/** The radix, or base of this system.  Set to 10000 */
 const RADIX: i32 = 10000;

/** The minimum exponent before underflow is signaled.  Flush to zero
     *  occurs at minExp-DIGITS */
 const MIN_EXP: i32 = -32767;

/** The maximum exponent before overflow is signaled and results flushed
     *  to infinity */
 const MAX_EXP: i32 = 32768;

/** The amount under/overflows are scaled by before going to trap handler */
 const ERR_SCALE: i32 = 32760;

/** Indicator value for normal finite numbers. */
 const FINITE: i8 = 0;

/** Indicator value for Infinity. */
 const INFINITE: i8 = 1;

/** Indicator value for signaling NaN. */
 const SNAN: i8 = 2;

/** Indicator value for quiet NaN. */
 const QNAN: i8 = 3;

/** String for NaN representation. */
 const NAN_STRING: &'static str = "NaN";

/** String for positive infinity representation. */
 const POS_INFINITY_STRING: &'static str = "Infinity";

/** String for negative infinity representation. */
 const NEG_INFINITY_STRING: &'static str = "-Infinity";

/** Name for traps triggered by addition. */
 const ADD_TRAP: &'static str = "add";

/** Name for traps triggered by multiplication. */
 const MULTIPLY_TRAP: &'static str = "multiply";

/** Name for traps triggered by division. */
 const DIVIDE_TRAP: &'static str = "divide";

/** Name for traps triggered by square root. */
 const SQRT_TRAP: &'static str = "sqrt";

/** Name for traps triggered by alignment. */
 const ALIGN_TRAP: &'static str = "align";

/** Name for traps triggered by truncation. */
 const TRUNC_TRAP: &'static str = "trunc";

/** Name for traps triggered by nextAfter. */
 const NEXT_AFTER_TRAP: &'static str = "nextAfter";

/** Name for traps triggered by lessThan. */
 const LESS_THAN_TRAP: &'static str = "lessThan";

/** Name for traps triggered by greaterThan. */
 const GREATER_THAN_TRAP: &'static str = "greaterThan";

/** Name for traps triggered by newInstance. */
 const NEW_INSTANCE_TRAP: &'static str = "newInstance";
 
pub struct Dfp {

    /** Mantissa. */
     pub mant: vec<i32>,
     
    /** Sign bit: 1 for positive, -1 for negative. */
     pub sign: i8,

    /** Exponent. */
     pub exp: i32,

    /** Indicator for non-finite / non-number values. */
     pub nans: i8,

    /** Factory building similar Dfp's. */
     pub field: DfpField
}

impl RealFieldElement<Dfp> for Dfp {
	
}

impl Dfp {

    /** Makes an instance with a value of zero.
     * @param field field to which this instance belongs
     */
    fn new( field: &DfpField) -> Dfp {
    	Dfp { 
        mant: std::Vec::new(),
        sign: 1,
        exp: 0,
        nans: FINITE,
        field: field }
    }

    /** Create an instance from a byte value.
     * @param field field to which this instance belongs
     * @param x value to convert to an instance
     */
    fn new_i8( field: &DfpField,  x: i8) -> Dfp {
        new_i64(field, x as i64);
    }

    /** Create an instance from an int value.
     * @param field field to which this instance belongs
     * @param x value to convert to an instance
     */
    fn new_i32( field: &DfpField,  x: i32) -> Dfp {
        new_i64(field, x as i64);
    }

    /** Create an instance from a long value.
     * @param field field to which this instance belongs
     * @param x value to convert to an instance
     */
    fn new_i64( field: &DfpField,  x: i64) -> Dfp {
        // initialize as if 0
        let mant = Vec::new();
        nans = FINITE;
        self.field = field;
        let is_long_min: bool = false;
        if x == i32::MIN {
            // special case for Long.MIN_VALUE (-9223372036854775808)
            // we must shift it before taking its absolute value
            is_long_min = true;
            x+= 1;
        }
        // set the sign
        if x < 0 {
            sign = -1;
            x = -x;
        } else {
            sign = 1;
        }
        exp = 0;
        while x != 0 {
            mant.push(x % RADIX) as i32;
            x /= RADIX;
            exp+= 1;
        }
        for i in exp..field.get_radix_digits() {
        	mant.insert(0,0);
        }
        if is_long_min {
            // we know in this case that fixing the last digit is sufficient
            for i in 0..mant.len() {
                if mant[i] != 0 {
                    mant[i]+= 1;
                    break;
                }
            }
        }
        Dfp { 
           	mant: mant,
        	sign: sign,
        	exp: exp,
        	nans: nans,
        	field: field 
        }
        
    }

    /** Create an instance from a double value.
     * @param field field to which this instance belongs
     * @param x value to convert to an instance
     */
    fn new_f64( field: &DfpField,  x: f64) -> Dfp {
        // initialize as if 0
       
         let bits: i64 = double_to_long_bits(&x);
         let mantissa: i64 = bits & 0x000fffffffffffff;
         let exponent: i32 = ((bits & 0x7ff0000000000000) >> 52) as i32 - 1023;
        if exponent == -1023 {
            // Zero or sub-normal
            if x == 0 {
                // make sure 0 has the right sign
                if (bits & 0x8000000000000000) != 0 {
                    sign = -1;
                }
                return;
            }
            exponent+= 1;
            // Normalize the subnormal number
            while (mantissa & 0x0010000000000000) == 0 {
                exponent-= 1;
                mantissa <<= 1;
            }
            mantissa &= 0x000fffffffffffff;
        }
        if exponent == 1024 {
            // infinity or NAN
            if x != x {
                sign = 1 as i8;
                nans = QNAN;
            } else if x < 0 {
                sign = -1 as i8;
                nans = INFINITE;
            } else {
                sign = 1 as i8;
                nans = INFINITE;
            }
            return;
        }
         let xdfp: Dfp = Dfp::new_i64(field, mantissa);
        xdfp = xdfp.divide(Dfp::new_i64(field, 4503599627370496)).add(field.get_one());
        xdfp = xdfp.multiply(DfpMath.pow(field.get_two(), exponent));
        if (bits & 0x8000000000000000) != 0 {
            xdfp = xdfp.negate();
        }
         Dfp { 
           	mant: xdfp.mant,
        	sign: xdfp.sign,
        	exp: xdfp.exp,
        	nans: dfp.nans,
        	field: field 
        }
        // Divide by 2^52, then add one
    }

    /** Copy constructor.
     * @param d instance to copy
     */
    pub fn new( d: &Dfp) -> Dfp {
        
        Dfp {
        	mant: d.mant,
        	sign: d.sign,
        	exp: d.exp,
        	nans: d.nans,
        	field: d.field 
        }
    }

    /** Create an instance from a String representation.
     * @param field field to which this instance belongs
     * @param s string representation of the instance
     */
    fn new_str( field: &DfpField,  s: &String) -> Dfp {
        // initialize as if 0
        let mant = Vec::new();
        let sign = 1;
        let exp = 0;
        let nans = FINITE;
         let decimal_found: bool = false;
        // size of radix in decimal digits
         let rsize: i32 = 4;
        // Starting offset into Striped
         let offset: i32 = 4;
         let striped: [Option<char>; get_radix_digits() * rsize + offset * 2] = [None; get_radix_digits() * rsize + offset * 2];
        // Check some special cases
        if s.equals(POS_INFINITY_STRING) {
            sign = 1 as i8;
            nans = INFINITE;
            return;
        } else if s.equals(NEG_INFINITY_STRING) {
            sign = -1 as i8;
            nans = INFINITE;
            return;
        } else if s.equals(NAN_STRING) {
            sign = 1 as i8;
            nans = QNAN;
            return;
        } else {
	        // Check for scientific notation
	         let p: i32 = s.index_of("e");
	        if p == -1 {
	            // try upper case?
	            p = s.index_of("E");
	        }
	         let fpdecimal: String;
	         let sciexp: i32 = 0;
	        if p != -1 {
	            // scientific notation
	            fpdecimal = s.substring(0, p);
	             let fpexp: String = s.substring(p + 1);
	             let negative: bool = false;
	             for i in 0..fpexp.len() {
	                if fpexp.char_at(i) == '-' {
	                    negative = true;
	                    continue;
	                }
	                if fpexp.char_at(i) >= '0' && fpexp.char_at(i) <= '9' {
	                    sciexp = sciexp * 10 + fpexp.char_at(i) - '0';
	                }
	            }
	            if negative {
	                sciexp = -sciexp;
	            }
	        } else {
	            // normal case
	            fpdecimal = s;
	        }
	        // If there is a minus sign in the number then it is negative
	        if fpdecimal.index_of("-") != -1 {
	            sign = -1;
	        }
	        // First off, find all of the leading zeros, trailing zeros, and significant digits
	        p = 0;
	        // Move p to first significant digit
	         let decimal_pos: i32 = 0;
	        loop {
	            if fpdecimal.char_at(p) >= '1' && fpdecimal.char_at(p) <= '9' {
	                break;
	            }
	            if decimal_found && fpdecimal.char_at(p) == '0' {
	                decimal_pos-= 1;
	            }
	            if fpdecimal.char_at(p) == '.' {
	                decimal_found = true;
	            }
	            p+= 1;
	            if p == fpdecimal.length() {
	                break;
	            }
	        }
	        // Copy the string onto Stripped
	         let q: i32 = offset;
	        striped[0] = '0';
	        striped[1] = '0';
	        striped[2] = '0';
	        striped[3] = '0';
	         let significant_digits: i32 = 0;
	        loop {
	            if p == (fpdecimal.length()) {
	                break;
	            }
	            // Don't want to run pass the end of the array
	            if q == mant.length * rsize + offset + 1 {
	                break;
	            }
	            if fpdecimal.char_at(p) == '.' {
	                decimal_found = true;
	                decimal_pos = significant_digits;
	                p+= 1;
	                continue;
	            }
	            if fpdecimal.char_at(p) < '0' || fpdecimal.char_at(p) > '9' {
	                p+= 1;
	                continue;
	            }
	            striped[q] = fpdecimal.char_at(p);
	            q+= 1;
	            p+= 1;
	            significant_digits+= 1;
	        }
	        // If the decimal point has been found then get rid of trailing zeros.
	        if decimal_found && q != offset {
	            loop {
	                q-= 1;
	                if q == offset {
	                    break;
	                }
	                if striped[q] == '0' {
	                    significant_digits-= 1;
	                } else {
	                    break;
	                }
	            }
	        }
	        // special case of numbers like "0.00000"
	        if decimal_found && significant_digits == 0 {
	            decimal_pos = 0;
	        }
	        // Implicit decimal point at end of number if not present
	        if !decimal_found {
	            decimal_pos = q - offset;
	        }
	        // Find the number of significant trailing zeros
	        // set q to point to first sig digit
	        q = offset;
	        p = significant_digits - 1 + offset;
	        while p > q {
	            if striped[p] != '0' {
	                break;
	            }
	            p-= 1;
	        }
	        // Make sure the decimal is on a mod 10000 boundary
	         let i: i32 = ((rsize * 100) - decimal_pos - sciexp % rsize) % rsize;
	        q -= i;
	        decimal_pos += i;
	        // Make the mantissa length right by adding zeros at the end if necessary
	        while (p - q) < (mant.length * rsize) {
	            for (i = 0; i < rsize; i+= 1) {
	                striped[p+= 1] = '0';
	            }
	        }
	        // and where the least significant digit is
	        for (i = mant.length - 1; i >= 0; i-= 1) {
	            mant[i] = (striped[q] - '0') * 1000 + (striped[q + 1] - '0') * 100 + (striped[q + 2] - '0') * 10 + (striped[q + 3] - '0');
	            q += 4;
	        }
	        exp = (decimal_pos + sciexp) / rsize;
	        if q < striped.length {
	            // Is there possible another digit?
	            round((striped[q] - '0') * 1000);
	        }
        }
    }

    /** Creates an instance with a non-finite value.
     * @param field field to which this instance belongs
     * @param sign sign of the Dfp to create
     * @param nans code of the value, must be one of {@link #INFINITE},
     * {@link #SNAN},  {@link #QNAN}
     */
    fn new( field: &DfpField,  sign: i8,  nans: i8) -> Dfp {
        self.field = field;
        self.mant = : [i32; field.get_radix_digits()] = [0; field.get_radix_digits()];
        self.sign = sign;
        self.exp = 0;
        self.nans = nans;
    }

    /** Create an instance with a value of 0.
     * Use this internally in preference to constructors to facilitate subclasses
     * @return a new instance with a value of 0
     */
    pub fn  new_instance(&self) -> Dfp  {
        return new Dfp(get_field());
    }

    /** Create an instance from a byte value.
     * @param x value to convert to an instance
     * @return a new instance with value x
     */
    pub fn  new_instance_i8(&self,  x: i8) -> Dfp  {
        return Dfp::new_i8(get_field(), x);
    }

    /** Create an instance from an int value.
     * @param x value to convert to an instance
     * @return a new instance with value x
     */
    pub fn  new_instance_i32(&self,  x: i32) -> Dfp  {
        return Dfp::new_i32(get_field(), x);
    }

    /** Create an instance from a long value.
     * @param x value to convert to an instance
     * @return a new instance with value x
     */
    pub fn  new_instance_i64(&self,  x: i64) -> Dfp  {
        return Dfp::new_i64(get_field(), x);
    }

    /** Create an instance from a double value.
     * @param x value to convert to an instance
     * @return a new instance with value x
     */
    pub fn  new_instance_f64(&self,  x: f64) -> Dfp  {
        return nDfp::new_f64(get_field(), x);
    }

    /** Create an instance by copying an existing one.
     * Use this internally in preference to constructors to facilitate subclasses.
     * @param d instance to copy
     * @return a new instance with the same value as d
     */
    pub fn  new_instance_dfp(&self,  d: &Dfp) -> Dfp  {
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != d.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            return dotrap(DfpField::FLAG_INVALID, NEW_INSTANCE_TRAP, d, result);
        }
        return new Dfp(d);
    }

    /** Create an instance from a String representation.
     * Use this internally in preference to constructors to facilitate subclasses.
     * @param s string representation of the instance
     * @return a new instance parsed from specified string
     */
    pub fn  new_instance_string(&self,  s: &String) -> Dfp  {
        return Dfp::new_string(field, s);
    }

    /** Creates an instance with a non-finite value.
     * @param sig sign of the Dfp to create
     * @param code code of the value, must be one of {@link #INFINITE},
     * {@link #SNAN},  {@link #QNAN}
     * @return a new instance with a non-finite value
     */
    pub fn  new_instance_sig_code(&self,  sig: i8,  code: i8) -> Dfp  {
        return field.new_dfp(sig, code);
    }

    /** Get the {@link org.apache.commons.math3.Field Field} (really a {@link DfpField}) to which the instance belongs.
     * <p>
     * The field is linked to the number of digits and acts as a factory
     * for {@link Dfp} instances.
     * </p>
     * @return {@link org.apache.commons.math3.Field Field} (really a {@link DfpField}) to which the instance belongs
     */
    pub fn  get_field(&self) -> DfpField  {
        return field;
    }

    /** Get the number of radix digits of the instance.
     * @return number of radix digits
     */
    pub fn  get_radix_digits(&self) -> i32  {
        return field.get_radix_digits();
    }

    /** Get the constant 0.
     * @return a Dfp with value zero
     */
    pub fn  get_zero(&self) -> Dfp  {
        return field.get_zero();
    }

    /** Get the constant 1.
     * @return a Dfp with value one
     */
    pub fn  get_one(&self) -> Dfp  {
        return field.get_one();
    }

    /** Get the constant 2.
     * @return a Dfp with value two
     */
    pub fn  get_two(&self) -> Dfp  {
        return field.get_two();
    }

    /** Shift the mantissa left, and adjust the exponent to compensate.
     */
    fn  shift_left(&self) -> void  {
        for ( let i: i32 = mant.length - 1; i > 0; i-= 1) {
            mant[i] = mant[i - 1];
        }
        mant[0] = 0;
        exp-= 1;
    }

    /* Note that shiftRight() does not call round() as that round() itself
     uses shiftRight() */
    /** Shift the mantissa right, and adjust the exponent to compensate.
     */
    fn  shift_right(&self) -> void  {
        for ( let i: i32 = 0; i < mant.length - 1; i+= 1) {
            mant[i] = mant[i + 1];
        }
        mant[mant.length - 1] = 0;
        exp+= 1;
    }

    /** Make our exp equal to the supplied one, this may cause rounding.
     *  Also causes de-normalized numbers.  These numbers are generally
     *  dangerous because most routines assume normalized numbers.
     *  Align doesn't round, so it will return the last digit destroyed
     *  by shifting right.
     *  @param e desired exponent
     *  @return last digit destroyed by shifting right
     */
    fn  align(&self,  e: i32) -> i32  {
         let lostdigit: i32 = 0;
         let inexact: bool = false;
         let diff: i32 = exp - e;
         let adiff: i32 = diff;
        if adiff < 0 {
            adiff = -adiff;
        }
        if diff == 0 {
            return 0;
        }
        if adiff > (mant.length + 1) {
            // Special case
            Arrays.fill(mant, 0);
            exp = e;
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
            dotrap(DfpField::FLAG_INEXACT, ALIGN_TRAP, self, self);
            return 0;
        }
        for ( let i: i32 = 0; i < adiff; i+= 1) {
            if diff < 0 {
                /* Keep track of loss -- only signal inexact after losing 2 digits.
                 * the first lost digit is returned to add() and may be incorporated
                 * into the result.
                 */
                if lostdigit != 0 {
                    inexact = true;
                }
                lostdigit = mant[0];
                shift_right();
            } else {
                shift_left();
            }
        }
        if inexact {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
            dotrap(DfpField::FLAG_INEXACT, ALIGN_TRAP, self, self);
        }
        return lostdigit;
    }

    /** Check if instance is less than x.
     * @param x number to check instance against
     * @return true if instance is less than x and neither are NaN, false otherwise
     */
    pub fn  less_than(&self,  x: &Dfp) -> bool  {
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != x.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, x, result);
            return false;
        }
        /* if a nan is involved, signal invalid and return false */
        if is_na_n() || x.is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, x, new_instance(get_zero()));
            return false;
        }
        return compare(self, x) < 0;
    }

    /** Check if instance is greater than x.
     * @param x number to check instance against
     * @return true if instance is greater than x and neither are NaN, false otherwise
     */
    pub fn  greater_than(&self,  x: &Dfp) -> bool  {
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != x.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            dotrap(DfpField::FLAG_INVALID, GREATER_THAN_TRAP, x, result);
            return false;
        }
        /* if a nan is involved, signal invalid and return false */
        if is_na_n() || x.is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, GREATER_THAN_TRAP, x, new_instance(get_zero()));
            return false;
        }
        return compare(self, x) > 0;
    }

    /** Check if instance is less than or equal to 0.
     * @return true if instance is not NaN and less than or equal to 0, false otherwise
     */
    pub fn  negative_or_null(&self) -> bool  {
        if is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, self, new_instance(get_zero()));
            return false;
        }
        return (sign < 0) || ((mant[mant.length - 1] == 0) && !is_infinite());
    }

    /** Check if instance is strictly less than 0.
     * @return true if instance is not NaN and less than or equal to 0, false otherwise
     */
    pub fn  strictly_negative(&self) -> bool  {
        if is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, self, new_instance(get_zero()));
            return false;
        }
        return (sign < 0) && ((mant[mant.length - 1] != 0) || is_infinite());
    }

    /** Check if instance is greater than or equal to 0.
     * @return true if instance is not NaN and greater than or equal to 0, false otherwise
     */
    pub fn  positive_or_null(&self) -> bool  {
        if is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, self, new_instance(get_zero()));
            return false;
        }
        return (sign > 0) || ((mant[mant.length - 1] == 0) && !is_infinite());
    }

    /** Check if instance is strictly greater than 0.
     * @return true if instance is not NaN and greater than or equal to 0, false otherwise
     */
    pub fn  strictly_positive(&self) -> bool  {
        if is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, self, new_instance(get_zero()));
            return false;
        }
        return (sign > 0) && ((mant[mant.length - 1] != 0) || is_infinite());
    }

    /** Get the absolute value of instance.
     * @return absolute value of instance
     * @since 3.2
     */
    pub fn  abs(&self) -> Dfp  {
         let result: Dfp = new_instance(self);
        result.sign = 1;
        return result;
    }

    /** Check if instance is infinite.
     * @return true if instance is infinite
     */
    pub fn  is_infinite(&self) -> bool  {
        return nans == INFINITE;
    }

    /** Check if instance is not a number.
     * @return true if instance is not a number
     */
    pub fn  is_nan(&self) -> bool  {
        return (nans == QNAN) || (nans == SNAN);
    }

    /** Check if instance is equal to zero.
     * @return true if instance is equal to zero
     */
    pub fn  is_zero(&self) -> bool  {
        if is_na_n() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            dotrap(DfpField::FLAG_INVALID, LESS_THAN_TRAP, self, new_instance(get_zero()));
            return false;
        }
        return (mant[mant.length - 1] == 0) && !is_infinite();
    }

    /** Check if instance is equal to x.
     * @param other object to check instance against
     * @return true if instance is equal to x and neither are NaN, false otherwise
     */
    pub fn  equals(&self,  other: &Object) -> bool  {
        if other instanceof Dfp {
             let x: Dfp = other as Dfp;
            if is_na_n() || x.is_na_n() || field.get_radix_digits() != x.field.get_radix_digits() {
                return false;
            }
            return compare(self, x) == 0;
        }
        return false;
    }

    /**
     * Gets a hashCode for the instance.
     * @return a hash code value for this object
     */
    pub fn  hash_code(&self) -> i32  {
        return 17 + ( if is_zero() { 0 } else { (sign << 8) }) + (nans << 16) + exp + Arrays.hash_code(mant);
    }

    /** Check if instance is not equal to x.
     * @param x number to check instance against
     * @return true if instance is not equal to x and neither are NaN, false otherwise
     */
    pub fn  unequal(&self,  x: &Dfp) -> bool  {
        if is_na_n() || x.is_na_n() || field.get_radix_digits() != x.field.get_radix_digits() {
            return false;
        }
        return greater_than(x) || less_than(x);
    }

    /** Compare two instances.
     * @param a first instance in comparison
     * @param b second instance in comparison
     * @return -1 if a<b, 1 if a>b and 0 if a==b
     *  Note this method does not properly handle NaNs or numbers with different precision.
     */
    fn  compare( a: &Dfp,  b: &Dfp) -> i32  {
        // Ignore the sign of zero
        if a.mant[a.mant.length - 1] == 0 && b.mant[b.mant.length - 1] == 0 && a.nans == FINITE && b.nans == FINITE {
            return 0;
        }
        if a.sign != b.sign {
            if a.sign == -1 {
                return -1;
            } else {
                return 1;
            }
        }
        // deal with the infinities
        if a.nans == INFINITE && b.nans == FINITE {
            return a.sign;
        }
        if a.nans == FINITE && b.nans == INFINITE {
            return -b.sign;
        }
        if a.nans == INFINITE && b.nans == INFINITE {
            return 0;
        }
        // Handle special case when a or b is zero, by ignoring the exponents
        if b.mant[b.mant.length - 1] != 0 && a.mant[b.mant.length - 1] != 0 {
            if a.exp < b.exp {
                return -a.sign;
            }
            if a.exp > b.exp {
                return a.sign;
            }
        }
        // compare the mantissas
        for ( let i: i32 = a.mant.length - 1; i >= 0; i-= 1) {
            if a.mant[i] > b.mant[i] {
                return a.sign;
            }
            if a.mant[i] < b.mant[i] {
                return -a.sign;
            }
        }
        return 0;
    }

    /** Round to nearest integer using the round-half-even method.
     *  That is round to nearest integer unless both are equidistant.
     *  In which case round to the even one.
     *  @return rounded value
     * @since 3.2
     */
    pub fn  rint(&self) -> Dfp  {
        return trunc(DfpField::RoundingMode::ROUND_HALF_EVEN);
    }

    /** Round to an integer using the round floor mode.
     * That is, round toward -Infinity
     *  @return rounded value
     * @since 3.2
     */
    pub fn  floor(&self) -> Dfp  {
        return trunc(DfpField::RoundingMode::ROUND_FLOOR);
    }

    /** Round to an integer using the round ceil mode.
     * That is, round toward +Infinity
     *  @return rounded value
     * @since 3.2
     */
    pub fn  ceil(&self) -> Dfp  {
        return trunc(DfpField::RoundingMode::ROUND_CEIL);
    }

    /** Returns the IEEE remainder.
     * @param d divisor
     * @return this less n × d, where n is the integer closest to this/d
     * @since 3.2
     */
    pub fn  remainder(&self,  d: &Dfp) -> Dfp  {
         let result: Dfp = self.subtract(self.divide(d).rint().multiply(d));
        // IEEE 854-1987 says that if the result is zero, then it carries the sign of this
        if result.mant[mant.length - 1] == 0 {
            result.sign = sign;
        }
        return result;
    }

    /** Does the integer conversions with the specified rounding.
     * @param rmode rounding mode to use
     * @return truncated value
     */
    fn  trunc(&self,  rmode: &DfpField.RoundingMode) -> Dfp  {
         let changed: bool = false;
        if is_na_n() {
            return new_instance(self);
        }
        if nans == INFINITE {
            return new_instance(self);
        }
        if mant[mant.length - 1] == 0 {
            // a is zero
            return new_instance(self);
        }
        /* If the exponent is less than zero then we can certainly
         * return zero */
        if exp < 0 {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
             let result: Dfp = new_instance(get_zero());
            result = dotrap(DfpField::FLAG_INEXACT, TRUNC_TRAP, self, result);
            return result;
        }
        if exp >= mant.length {
            return new_instance(self);
        }
        /* General case:  create another dfp, result, that contains the
         * a with the fractional part lopped off.  */
         let result: Dfp = new_instance(self);
        for ( let i: i32 = 0; i < mant.length - result.exp; i+= 1) {
            changed |= result.mant[i] != 0;
            result.mant[i] = 0;
        }
        if changed {
            switch(rmode) {
                case ROUND_FLOOR:
                    if result.sign == -1 {
                        // then we must increment the mantissa by one
                        result = result.add(new_instance(-1));
                    }
                    break;
                case ROUND_CEIL:
                    if result.sign == 1 {
                        // then we must increment the mantissa by one
                        result = result.add(get_one());
                    }
                    break;
                case ROUND_HALF_EVEN:
                default:
                     let half: Dfp = new_instance("0.5");
                    // difference between this and result
                     let a: Dfp = subtract(result);
                    // force positive (take abs)
                    a.sign = 1;
                    if a.greater_than(half) {
                        a = new_instance(get_one());
                        a.sign = sign;
                        result = result.add(a);
                    }
                    /** If exactly equal to 1/2 and odd then increment */
                    if a.equals(half) && result.exp > 0 && (result.mant[mant.length - result.exp] & 1) != 0 {
                        a = new_instance(get_one());
                        a.sign = sign;
                        result = result.add(a);
                    }
                    break;
            }
            // signal inexact
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
            result = dotrap(DfpField::FLAG_INEXACT, TRUNC_TRAP, self, result);
            return result;
        }
        return result;
    }

    /** Convert this to an integer.
     * If greater than 2147483647, it returns 2147483647. If less than -2147483648 it returns -2147483648.
     * @return converted number
     */
    pub fn  int_value(&self) -> i32  {
         let rounded: Dfp;
         let result: i32 = 0;
        rounded = rint();
        if rounded.greater_than(new_instance(2147483647)) {
            return 2147483647;
        }
        if rounded.less_than(new_instance(-2147483648)) {
            return -2147483648;
        }
        for ( let i: i32 = mant.length - 1; i >= mant.length - rounded.exp; i-= 1) {
            result = result * RADIX + rounded.mant[i];
        }
        if rounded.sign == -1 {
            result = -result;
        }
        return result;
    }

    /** Get the exponent of the greatest power of 10000 that is
     *  less than or equal to the absolute value of this.  I.E.  if
     *  this is 10<sup>6</sup> then log10K would return 1.
     *  @return integer base 10000 logarithm
     */
    pub fn  log10_k(&self) -> i32  {
        return exp - 1;
    }

    /** Get the specified  power of 10000.
     * @param e desired power
     * @return 10000<sup>e</sup>
     */
    pub fn  power10_k(&self,  e: i32) -> Dfp  {
         let d: Dfp = new_instance(get_one());
        d.exp = e + 1;
        return d;
    }

    /** Get the exponent of the greatest power of 10 that is less than or equal to abs(this).
     *  @return integer base 10 logarithm
     * @since 3.2
     */
    pub fn  int_log10(&self) -> i32  {
        if mant[mant.length - 1] > 1000 {
            return exp * 4 - 1;
        }
        if mant[mant.length - 1] > 100 {
            return exp * 4 - 2;
        }
        if mant[mant.length - 1] > 10 {
            return exp * 4 - 3;
        }
        return exp * 4 - 4;
    }

    /** Return the specified  power of 10.
     * @param e desired power
     * @return 10<sup>e</sup>
     */
    pub fn  power10(&self,  e: i32) -> Dfp  {
         let d: Dfp = new_instance(get_one());
        if e >= 0 {
            d.exp = e / 4 + 1;
        } else {
            d.exp = (e + 1) / 4;
        }
        switch((e % 4 + 4) % 4) {
            case 0:
                break;
            case 1:
                d = d.multiply(10);
                break;
            case 2:
                d = d.multiply(100);
                break;
            default:
                d = d.multiply(1000);
        }
        return d;
    }

    /** Negate the mantissa of this by computing the complement.
     *  Leaves the sign bit unchanged, used internally by add.
     *  Denormalized numbers are handled properly here.
     *  @param extra ???
     *  @return ???
     */
    fn  complement(&self,  extra: i32) -> i32  {
        extra = RADIX - extra;
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
            mant[i] = RADIX - mant[i] - 1;
        }
         let rh: i32 = extra / RADIX;
        extra -= rh * RADIX;
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
             let r: i32 = mant[i] + rh;
            rh = r / RADIX;
            mant[i] = r - rh * RADIX;
        }
        return extra;
    }

    /** Add x to this.
     * @param x number to add
     * @return sum of this and x
     */
    pub fn  add(&self,  x: &Dfp) -> Dfp  {
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != x.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            return dotrap(DfpField::FLAG_INVALID, ADD_TRAP, x, result);
        }
        /* handle special cases */
        if nans != FINITE || x.nans != FINITE {
            if is_na_n() {
                return self;
            }
            if x.is_na_n() {
                return x;
            }
            if nans == INFINITE && x.nans == FINITE {
                return self;
            }
            if x.nans == INFINITE && nans == FINITE {
                return x;
            }
            if x.nans == INFINITE && nans == INFINITE && sign == x.sign {
                return x;
            }
            if x.nans == INFINITE && nans == INFINITE && sign != x.sign {
                field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
                 let result: Dfp = new_instance(get_zero());
                result.nans = QNAN;
                result = dotrap(DfpField::FLAG_INVALID, ADD_TRAP, x, result);
                return result;
            }
        }
        /* copy this and the arg */
         let a: Dfp = new_instance(self);
         let b: Dfp = new_instance(x);
        /* initialize the result object */
         let result: Dfp = new_instance(get_zero());
        /* Make all numbers positive, but remember their sign */
         let asign: i8 = a.sign;
         let bsign: i8 = b.sign;
        a.sign = 1;
        b.sign = 1;
        /* The result will be signed like the arg with greatest magnitude */
         let rsign: i8 = bsign;
        if compare(a, b) > 0 {
            rsign = asign;
        }
        /* Handle special case when a or b is zero, by setting the exponent
       of the zero number equal to the other one.  This avoids an alignment
       which would cause catastropic loss of precision */
        if b.mant[mant.length - 1] == 0 {
            b.exp = a.exp;
        }
        if a.mant[mant.length - 1] == 0 {
            a.exp = b.exp;
        }
        /* align number with the smaller exponent */
         let aextradigit: i32 = 0;
         let bextradigit: i32 = 0;
        if a.exp < b.exp {
            aextradigit = a.align(b.exp);
        } else {
            bextradigit = b.align(a.exp);
        }
        /* complement the smaller of the two if the signs are different */
        if asign != bsign {
            if asign == rsign {
                bextradigit = b.complement(bextradigit);
            } else {
                aextradigit = a.complement(aextradigit);
            }
        }
        /* add the mantissas */
         let rh: i32 = 0;
        /* acts as a carry */
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
             let r: i32 = a.mant[i] + b.mant[i] + rh;
            rh = r / RADIX;
            result.mant[i] = r - rh * RADIX;
        }
        result.exp = a.exp;
        result.sign = rsign;
        if rh != 0 && (asign == bsign) {
             let lostdigit: i32 = result.mant[0];
            result.shift_right();
            result.mant[mant.length - 1] = rh;
             let excp: i32 = result.round(lostdigit);
            if excp != 0 {
                result = dotrap(excp, ADD_TRAP, x, result);
            }
        }
        /* normalize the result */
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
            if result.mant[mant.length - 1] != 0 {
                break;
            }
            result.shift_left();
            if i == 0 {
                result.mant[0] = aextradigit + bextradigit;
                aextradigit = 0;
                bextradigit = 0;
            }
        }
        /* result is zero if after normalization the most sig. digit is zero */
        if result.mant[mant.length - 1] == 0 {
            result.exp = 0;
            if asign != bsign {
                // Unless adding 2 negative zeros, sign is positive
                // Per IEEE 854-1987 Section 6.3
                result.sign = 1;
            }
        }
        /* Call round to test for over/under flows */
         let excp: i32 = result.round(aextradigit + bextradigit);
        if excp != 0 {
            result = dotrap(excp, ADD_TRAP, x, result);
        }
        return result;
    }

    /** Returns a number that is this number with the sign bit reversed.
     * @return the opposite of this
     */
    pub fn  negate(&self) -> Dfp  {
         let result: Dfp = new_instance(self);
        result.sign = -result.sign as i8;
        return result;
    }

    /** Subtract x from this.
     * @param x number to subtract
     * @return difference of this and a
     */
    pub fn  subtract(&self,  x: &Dfp) -> Dfp  {
        return add(x.negate());
    }

    /** Round this given the next digit n using the current rounding mode.
     * @param n ???
     * @return the IEEE flag if an exception occurred
     */
    fn  round(&self,  n: i32) -> i32  {
         let inc: bool = false;
        switch(field.get_rounding_mode()) {
            case ROUND_DOWN:
                inc = false;
                break;
            case ROUND_UP:
                // round up if n!=0
                inc = n != 0;
                break;
            case ROUND_HALF_UP:
                // round half up
                inc = n >= 5000;
                break;
            case ROUND_HALF_DOWN:
                // round half down
                inc = n > 5000;
                break;
            case ROUND_HALF_EVEN:
                // round half-even
                inc = n > 5000 || (n == 5000 && (mant[0] & 1) == 1);
                break;
            case ROUND_HALF_ODD:
                // round half-odd
                inc = n > 5000 || (n == 5000 && (mant[0] & 1) == 0);
                break;
            case ROUND_CEIL:
                // round ceil
                inc = sign == 1 && n != 0;
                break;
            case ROUND_FLOOR:
            default:
                // round floor
                inc = sign == -1 && n != 0;
                break;
        }
        if inc {
            // increment if necessary
             let rh: i32 = 1;
            for ( let i: i32 = 0; i < mant.length; i+= 1) {
                 let r: i32 = mant[i] + rh;
                rh = r / RADIX;
                mant[i] = r - rh * RADIX;
            }
            if rh != 0 {
                shift_right();
                mant[mant.length - 1] = rh;
            }
        }
        // check for exceptional cases and raise signals if necessary
        if exp < MIN_EXP {
            // Gradual Underflow
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_UNDERFLOW);
            return DfpField::FLAG_UNDERFLOW;
        }
        if exp > MAX_EXP {
            // Overflow
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_OVERFLOW);
            return DfpField::FLAG_OVERFLOW;
        }
        if n != 0 {
            // Inexact
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
            return DfpField::FLAG_INEXACT;
        }
        return 0;
    }

    /** Multiply this by x.
     * @param x multiplicand
     * @return product of this and x
     */
    pub fn  multiply(&self,  x: &Dfp) -> Dfp  {
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != x.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            return dotrap(DfpField::FLAG_INVALID, MULTIPLY_TRAP, x, result);
        }
         let result: Dfp = new_instance(get_zero());
        /* handle special cases */
        if nans != FINITE || x.nans != FINITE {
            if is_na_n() {
                return self;
            }
            if x.is_na_n() {
                return x;
            }
            if nans == INFINITE && x.nans == FINITE && x.mant[mant.length - 1] != 0 {
                result = new_instance(self);
                result.sign = (sign * x.sign) as i8;
                return result;
            }
            if x.nans == INFINITE && nans == FINITE && mant[mant.length - 1] != 0 {
                result = new_instance(x);
                result.sign = (sign * x.sign) as i8;
                return result;
            }
            if x.nans == INFINITE && nans == INFINITE {
                result = new_instance(self);
                result.sign = (sign * x.sign) as i8;
                return result;
            }
            if (x.nans == INFINITE && nans == FINITE && mant[mant.length - 1] == 0) || (nans == INFINITE && x.nans == FINITE && x.mant[mant.length - 1] == 0) {
                field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
                result = new_instance(get_zero());
                result.nans = QNAN;
                result = dotrap(DfpField::FLAG_INVALID, MULTIPLY_TRAP, x, result);
                return result;
            }
        }
        // Big enough to hold even the largest result
         let product: [i32; mant.length * 2] = [0; mant.length * 2];
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
            // acts as a carry
             let rh: i32 = 0;
            for ( let j: i32 = 0; j < mant.length; j+= 1) {
                // multiply the 2 digits
                 let r: i32 = mant[i] * x.mant[j];
                // add to the product digit with carry in
                r += product[i + j] + rh;
                rh = r / RADIX;
                product[i + j] = r - rh * RADIX;
            }
            product[i + mant.length] = rh;
        }
        // Find the most sig digit
        // default, in case result is zero
         let md: i32 = mant.length * 2 - 1;
        for ( let i: i32 = mant.length * 2 - 1; i >= 0; i-= 1) {
            if product[i] != 0 {
                md = i;
                break;
            }
        }
        // Copy the digits into the result
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
            result.mant[mant.length - i - 1] = product[md - i];
        }
        // Fixup the exponent.
        result.exp = exp + x.exp + md - 2 * mant.length + 1;
        result.sign = ( if (sign == x.sign) { 1 } else { -1 }) as i8;
        if result.mant[mant.length - 1] == 0 {
            // if result is zero, set exp to zero
            result.exp = 0;
        }
         let excp: i32;
        if md > (mant.length - 1) {
            excp = result.round(product[md - mant.length]);
        } else {
            // has no effect except to check status
            excp = result.round(0);
        }
        if excp != 0 {
            result = dotrap(excp, MULTIPLY_TRAP, x, result);
        }
        return result;
    }

    /** Multiply this by a single digit x.
     * @param x multiplicand
     * @return product of this and x
     */
    pub fn  multiply(&self,  x: i32) -> Dfp  {
        if x >= 0 && x < RADIX {
            return multiply_fast(x);
        } else {
            return multiply(new_instance(x));
        }
    }

    /** Multiply this by a single digit 0<=x<radix.
     * There are speed advantages in this special case.
     * @param x multiplicand
     * @return product of this and x
     */
    fn  multiply_fast(&self,  x: i32) -> Dfp  {
         let result: Dfp = new_instance(self);
        /* handle special cases */
        if nans != FINITE {
            if is_na_n() {
                return self;
            }
            if nans == INFINITE && x != 0 {
                result = new_instance(self);
                return result;
            }
            if nans == INFINITE && x == 0 {
                field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
                result = new_instance(get_zero());
                result.nans = QNAN;
                result = dotrap(DfpField::FLAG_INVALID, MULTIPLY_TRAP, new_instance(get_zero()), result);
                return result;
            }
        }
        /* range check x */
        if x < 0 || x >= RADIX {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            result = new_instance(get_zero());
            result.nans = QNAN;
            result = dotrap(DfpField::FLAG_INVALID, MULTIPLY_TRAP, result, result);
            return result;
        }
         let rh: i32 = 0;
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
             let r: i32 = mant[i] * x + rh;
            rh = r / RADIX;
            result.mant[i] = r - rh * RADIX;
        }
         let lostdigit: i32 = 0;
        if rh != 0 {
            lostdigit = result.mant[0];
            result.shift_right();
            result.mant[mant.length - 1] = rh;
        }
        if result.mant[mant.length - 1] == 0 {
            // if result is zero, set exp to zero
            result.exp = 0;
        }
         let excp: i32 = result.round(lostdigit);
        if excp != 0 {
            result = dotrap(excp, MULTIPLY_TRAP, result, result);
        }
        return result;
    }

    /** Divide this by divisor.
     * @param divisor divisor
     * @return quotient of this by divisor
     */
    pub fn  divide(&self,  divisor: &Dfp) -> Dfp  {
        // current status of the dividend
         let dividend: i32;
        // quotient
         let quotient: i32;
        // remainder
         let remainder: i32;
        // current quotient digit we're working with
         let qd: i32;
        // number of significant quotient digits we have
         let nsqd: i32;
        // trial quotient digit
         let trial: i32 = 0;
        // minimum adjustment
         let minadj: i32;
        // Flag to indicate a good trail digit
         let trialgood: bool;
        // most sig digit in result
         let md: i32 = 0;
        // exceptions
         let excp: i32;
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != divisor.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            return dotrap(DfpField::FLAG_INVALID, DIVIDE_TRAP, divisor, result);
        }
         let result: Dfp = new_instance(get_zero());
        /* handle special cases */
        if nans != FINITE || divisor.nans != FINITE {
            if is_na_n() {
                return self;
            }
            if divisor.is_na_n() {
                return divisor;
            }
            if nans == INFINITE && divisor.nans == FINITE {
                result = new_instance(self);
                result.sign = (sign * divisor.sign) as i8;
                return result;
            }
            if divisor.nans == INFINITE && nans == FINITE {
                result = new_instance(get_zero());
                result.sign = (sign * divisor.sign) as i8;
                return result;
            }
            if divisor.nans == INFINITE && nans == INFINITE {
                field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
                result = new_instance(get_zero());
                result.nans = QNAN;
                result = dotrap(DfpField::FLAG_INVALID, DIVIDE_TRAP, divisor, result);
                return result;
            }
        }
        /* Test for divide by zero */
        if divisor.mant[mant.length - 1] == 0 {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_DIV_ZERO);
            result = new_instance(get_zero());
            result.sign = (sign * divisor.sign) as i8;
            result.nans = INFINITE;
            result = dotrap(DfpField::FLAG_DIV_ZERO, DIVIDE_TRAP, divisor, result);
            return result;
        }
        // one extra digit needed
        dividend = : [i32; mant.length + 1] = [0; mant.length + 1];
        // two extra digits needed 1 for overflow, 1 for rounding
        quotient = : [i32; mant.length + 2] = [0; mant.length + 2];
        // one extra digit needed
        remainder = : [i32; mant.length + 1] = [0; mant.length + 1];
        /* Initialize our most significant digits to zero */
        dividend[mant.length] = 0;
        quotient[mant.length] = 0;
        quotient[mant.length + 1] = 0;
        remainder[mant.length] = 0;
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
            dividend[i] = mant[i];
            quotient[i] = 0;
            remainder[i] = 0;
        }
        /* outer loop.  Once per quotient digit */
        nsqd = 0;
        for (qd = mant.length + 1; qd >= 0; qd-= 1) {
            /* Determine outer limits of our quotient digit */
            // r =  most sig 2 digits of dividend
             let div_msb: i32 = dividend[mant.length] * RADIX + dividend[mant.length - 1];
             let min: i32 = div_msb / (divisor.mant[mant.length - 1] + 1);
             let max: i32 = (div_msb + 1) / divisor.mant[mant.length - 1];
            trialgood = false;
            while !trialgood {
                // try the mean
                trial = (min + max) / 2;
                /* Multiply by divisor and store as remainder */
                 let rh: i32 = 0;
                for ( let i: i32 = 0; i < mant.length + 1; i+= 1) {
                     let dm: i32 =  if (i < mant.length) { divisor.mant[i] } else { 0 };
                     let r: i32 = (dm * trial) + rh;
                    rh = r / RADIX;
                    remainder[i] = r - rh * RADIX;
                }
                /* subtract the remainder from the dividend */
                // carry in to aid the subtraction
                rh = 1;
                for ( let i: i32 = 0; i < mant.length + 1; i+= 1) {
                     let r: i32 = ((RADIX - 1) - remainder[i]) + dividend[i] + rh;
                    rh = r / RADIX;
                    remainder[i] = r - rh * RADIX;
                }
                /* Lets analyze what we have here */
                if rh == 0 {
                    // trial is too big -- negative remainder
                    max = trial - 1;
                    continue;
                }
                /* find out how far off the remainder is telling us we are */
                minadj = (remainder[mant.length] * RADIX) + remainder[mant.length - 1];
                minadj /= divisor.mant[mant.length - 1] + 1;
                if minadj >= 2 {
                    // update the minimum
                    min = trial + minadj;
                    continue;
                }
                /* May have a good one here, check more thoroughly.  Basically
           its a good one if it is less than the divisor */
                // assume false
                trialgood = false;
                for ( let i: i32 = mant.length - 1; i >= 0; i-= 1) {
                    if divisor.mant[i] > remainder[i] {
                        trialgood = true;
                    }
                    if divisor.mant[i] < remainder[i] {
                        break;
                    }
                }
                if remainder[mant.length] != 0 {
                    trialgood = false;
                }
                if trialgood == false {
                    min = trial + 1;
                }
            }
            /* Great we have a digit! */
            quotient[qd] = trial;
            if trial != 0 || nsqd != 0 {
                nsqd+= 1;
            }
            if field.get_rounding_mode() == DfpField::RoundingMode::ROUND_DOWN && nsqd == mant.length {
                // We have enough for this mode
                break;
            }
            if nsqd > mant.length {
                // We have enough digits
                break;
            }
            /* move the remainder into the dividend while left shifting */
            dividend[0] = 0;
            for ( let i: i32 = 0; i < mant.length; i+= 1) {
                dividend[i + 1] = remainder[i];
            }
        }
        /* Find the most sig digit */
        // default
        md = mant.length;
        for ( let i: i32 = mant.length + 1; i >= 0; i-= 1) {
            if quotient[i] != 0 {
                md = i;
                break;
            }
        }
        /* Copy the digits into the result */
        for ( let i: i32 = 0; i < mant.length; i+= 1) {
            result.mant[mant.length - i - 1] = quotient[md - i];
        }
        /* Fixup the exponent. */
        result.exp = exp - divisor.exp + md - mant.length;
        result.sign = ( if (sign == divisor.sign) { 1 } else { -1 }) as i8;
        if result.mant[mant.length - 1] == 0 {
            // if result is zero, set exp to zero
            result.exp = 0;
        }
        if md > (mant.length - 1) {
            excp = result.round(quotient[md - mant.length]);
        } else {
            excp = result.round(0);
        }
        if excp != 0 {
            result = dotrap(excp, DIVIDE_TRAP, divisor, result);
        }
        return result;
    }

    /** Divide by a single digit less than radix.
     *  Special case, so there are speed advantages. 0 <= divisor < radix
     * @param divisor divisor
     * @return quotient of this by divisor
     */
    pub fn  divide(&self,  divisor: i32) -> Dfp  {
        // Handle special cases
        if nans != FINITE {
            if is_na_n() {
                return self;
            }
            if nans == INFINITE {
                return new_instance(self);
            }
        }
        // Test for divide by zero
        if divisor == 0 {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_DIV_ZERO);
             let result: Dfp = new_instance(get_zero());
            result.sign = sign;
            result.nans = INFINITE;
            result = dotrap(DfpField::FLAG_DIV_ZERO, DIVIDE_TRAP, get_zero(), result);
            return result;
        }
        // range check divisor
        if divisor < 0 || divisor >= RADIX {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            result = dotrap(DfpField::FLAG_INVALID, DIVIDE_TRAP, result, result);
            return result;
        }
         let result: Dfp = new_instance(self);
         let rl: i32 = 0;
        for ( let i: i32 = mant.length - 1; i >= 0; i-= 1) {
             let r: i32 = rl * RADIX + result.mant[i];
             let rh: i32 = r / divisor;
            rl = r - rh * divisor;
            result.mant[i] = rh;
        }
        if result.mant[mant.length - 1] == 0 {
            // normalize
            result.shift_left();
            // compute the next digit and put it in
             let r: i32 = rl * RADIX;
             let rh: i32 = r / divisor;
            rl = r - rh * divisor;
            result.mant[0] = rh;
        }
        // do the rounding
         let excp: i32 = result.round(rl * RADIX / divisor);
        if excp != 0 {
            result = dotrap(excp, DIVIDE_TRAP, result, result);
        }
        return result;
    }

    /** {@inheritDoc} */
    pub fn  reciprocal(&self) -> Dfp  {
        return field.get_one().divide(self);
    }

    /** Compute the square root.
     * @return square root of the instance
     * @since 3.2
     */
    pub fn  sqrt(&self) -> Dfp  {
        // check for unusual cases
        if nans == FINITE && mant[mant.length - 1] == 0 {
            // if zero
            return new_instance(self);
        }
        if nans != FINITE {
            if nans == INFINITE && sign == 1 {
                // if positive infinity
                return new_instance(self);
            }
            if nans == QNAN {
                return new_instance(self);
            }
            if nans == SNAN {
                 let result: Dfp;
                field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
                result = new_instance(self);
                result = dotrap(DfpField::FLAG_INVALID, SQRT_TRAP, null, result);
                return result;
            }
        }
        if sign == -1 {
            // if negative
             let result: Dfp;
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
            result = new_instance(self);
            result.nans = QNAN;
            result = dotrap(DfpField::FLAG_INVALID, SQRT_TRAP, null, result);
            return result;
        }
         let x: Dfp = new_instance(self);
        /* Lets make a reasonable guess as to the size of the square root */
        if x.exp < -1 || x.exp > 1 {
            x.exp = self.exp / 2;
        }
        /* Coarsely estimate the mantissa */
        switch(x.mant[mant.length - 1] / 2000) {
            case 0:
                x.mant[mant.length - 1] = x.mant[mant.length - 1] / 2 + 1;
                break;
            case 2:
                x.mant[mant.length - 1] = 1500;
                break;
            case 3:
                x.mant[mant.length - 1] = 2200;
                break;
            default:
                x.mant[mant.length - 1] = 3000;
        }
         let dx: Dfp = new_instance(x);
        /* Now that we have the first pass estimate, compute the rest
       by the formula dx = (y - x*x) / (2x); */
         let px: Dfp = get_zero();
         let ppx: Dfp = get_zero();
        while x.unequal(px) {
            dx = new_instance(x);
            dx.sign = -1;
            dx = dx.add(self.divide(x));
            dx = dx.divide(2);
            ppx = px;
            px = x;
            x = x.add(dx);
            if x.equals(ppx) {
                // alternating between two values
                break;
            }
            // is a sufficient test since dx is normalized
            if dx.mant[mant.length - 1] == 0 {
                break;
            }
        }
        return x;
    }

    /** Get a string representation of the instance.
     * @return string representation of the instance
     */
    pub fn  to_string(&self) -> String  {
        if nans != FINITE {
            // if non-finite exceptional cases
            if nans == INFINITE {
                return  if (sign < 0) { NEG_INFINITY_STRING } else { POS_INFINITY_STRING };
            } else {
                return NAN_STRING;
            }
        }
        if exp > mant.length || exp < -1 {
            return dfp2sci();
        }
        return dfp2string();
    }

    /** Convert an instance to a string using scientific notation.
     * @return string representation of the instance in scientific notation
     */
    fn  dfp2sci(&self) -> String  {
         let rawdigits: [Option<char>; mant.length * 4] = [None; mant.length * 4];
         let outputbuffer: [Option<char>; mant.length * 4 + 20] = [None; mant.length * 4 + 20];
         let p: i32;
         let q: i32;
         let e: i32;
         let ae: i32;
         let shf: i32;
        // Get all the digits
        p = 0;
        for ( let i: i32 = mant.length - 1; i >= 0; i-= 1) {
            rawdigits[p+= 1] = ((mant[i] / 1000) + '0') as char;
            rawdigits[p+= 1] = (((mant[i] / 100) % 10) + '0') as char;
            rawdigits[p+= 1] = (((mant[i] / 10) % 10) + '0') as char;
            rawdigits[p+= 1] = (((mant[i]) % 10) + '0') as char;
        }
        // Find the first non-zero one
        for (p = 0; p < rawdigits.length; p+= 1) {
            if rawdigits[p] != '0' {
                break;
            }
        }
        shf = p;
        // Now do the conversion
        q = 0;
        if sign == -1 {
            outputbuffer[q+= 1] = '-';
        }
        if p != rawdigits.length {
            // there are non zero digits...
            outputbuffer[q+= 1] = rawdigits[p+= 1];
            outputbuffer[q+= 1] = '.';
            while p < rawdigits.length {
                outputbuffer[q+= 1] = rawdigits[p+= 1];
            }
        } else {
            outputbuffer[q+= 1] = '0';
            outputbuffer[q+= 1] = '.';
            outputbuffer[q+= 1] = '0';
            outputbuffer[q+= 1] = 'e';
            outputbuffer[q+= 1] = '0';
            return new String(outputbuffer, 0, 5);
        }
        outputbuffer[q+= 1] = 'e';
        // Find the msd of the exponent
        e = exp * 4 - shf - 1;
        ae = e;
        if e < 0 {
            ae = -e;
        }
        // Find the largest p such that p < e
        for (p = 1000000000; p > ae; p /= 10) {
        // nothing to do
        }
        if e < 0 {
            outputbuffer[q+= 1] = '-';
        }
        while p > 0 {
            outputbuffer[q+= 1] = (ae / p + '0') as char;
            ae %= p;
            p /= 10;
        }
        return new String(outputbuffer, 0, q);
    }

    /** Convert an instance to a string using normal notation.
     * @return string representation of the instance in normal notation
     */
    fn  dfp2string(&self) -> String  {
         let buffer: [Option<char>; mant.length * 4 + 20] = [None; mant.length * 4 + 20];
         let p: i32 = 1;
         let q: i32;
         let e: i32 = exp;
         let point_inserted: bool = false;
        buffer[0] = ' ';
        if e <= 0 {
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '.';
            point_inserted = true;
        }
        while e < 0 {
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '0';
            e+= 1;
        }
        for ( let i: i32 = mant.length - 1; i >= 0; i-= 1) {
            buffer[p+= 1] = ((mant[i] / 1000) + '0') as char;
            buffer[p+= 1] = (((mant[i] / 100) % 10) + '0') as char;
            buffer[p+= 1] = (((mant[i] / 10) % 10) + '0') as char;
            buffer[p+= 1] = (((mant[i]) % 10) + '0') as char;
            if e-= 1 == 0 {
                buffer[p+= 1] = '.';
                point_inserted = true;
            }
        }
        while e > 0 {
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '0';
            buffer[p+= 1] = '0';
            e-= 1;
        }
        if !point_inserted {
            // Ensure we have a radix point!
            buffer[p+= 1] = '.';
        }
        // Suppress leading zeros
        q = 1;
        while buffer[q] == '0' {
            q+= 1;
        }
        if buffer[q] == '.' {
            q-= 1;
        }
        // Suppress trailing zeros
        while buffer[p - 1] == '0' {
            p-= 1;
        }
        // Insert sign
        if sign < 0 {
            buffer[q-= 1] = '-';
        }
        return new String(buffer, q, p - q);
    }

    /** Raises a trap.  This does not set the corresponding flag however.
     *  @param type the trap type
     *  @param what - name of routine trap occurred in
     *  @param oper - input operator to function
     *  @param result - the result computed prior to the trap
     *  @return The suggested return value from the trap handler
     */
    pub fn  dotrap(&self,  type: i32,  what: &String,  oper: &Dfp,  result: &Dfp) -> Dfp  {
         let def: Dfp = result;
        switch(type) {
            case DfpField::FLAG_INVALID:
                def = new_instance(get_zero());
                def.sign = result.sign;
                def.nans = QNAN;
                break;
            case DfpField::FLAG_DIV_ZERO:
                if nans == FINITE && mant[mant.length - 1] != 0 {
                    // normal case, we are finite, non-zero
                    def = new_instance(get_zero());
                    def.sign = (sign * oper.sign) as i8;
                    def.nans = INFINITE;
                }
                if nans == FINITE && mant[mant.length - 1] == 0 {
                    //  0/0
                    def = new_instance(get_zero());
                    def.nans = QNAN;
                }
                if nans == INFINITE || nans == QNAN {
                    def = new_instance(get_zero());
                    def.nans = QNAN;
                }
                if nans == INFINITE || nans == SNAN {
                    def = new_instance(get_zero());
                    def.nans = QNAN;
                }
                break;
            case DfpField::FLAG_UNDERFLOW:
                if (result.exp + mant.length) < MIN_EXP {
                    def = new_instance(get_zero());
                    def.sign = result.sign;
                } else {
                    // gradual underflow
                    def = new_instance(result);
                }
                result.exp += ERR_SCALE;
                break;
            case DfpField::FLAG_OVERFLOW:
                result.exp -= ERR_SCALE;
                def = new_instance(get_zero());
                def.sign = result.sign;
                def.nans = INFINITE;
                break;
            default:
                def = result;
                break;
        }
        return trap(type, what, oper, def, result);
    }

    /** Trap handler.  Subclasses may override this to provide trap
     *  functionality per IEEE 854-1987.
     *
     *  @param type  The exception type - e.g. FLAG_OVERFLOW
     *  @param what  The name of the routine we were in e.g. divide()
     *  @param oper  An operand to this function if any
     *  @param def   The default return value if trap not enabled
     *  @param result    The result that is specified to be delivered per
     *                   IEEE 854, if any
     *  @return the value that should be return by the operation triggering the trap
     */
    fn  trap(&self,  type: i32,  what: &String,  oper: &Dfp,  def: &Dfp,  result: &Dfp) -> Dfp  {
        return def;
    }

    /** Returns the type - one of FINITE, INFINITE, SNAN, QNAN.
     * @return type of the number
     */
    pub fn  classify(&self) -> i32  {
        return nans;
    }

    /** Creates an instance that is the same as x except that it has the sign of y.
     * abs(x) = dfp.copysign(x, dfp.one)
     * @param x number to get the value from
     * @param y number to get the sign from
     * @return a number with the value of x and the sign of y
     */
    pub fn  copysign( x: &Dfp,  y: &Dfp) -> Dfp  {
         let result: Dfp = x.new_instance(x);
        result.sign = y.sign;
        return result;
    }

    /** Returns the next number greater than this one in the direction of x.
     * If this==x then simply returns this.
     * @param x direction where to look at
     * @return closest number next to instance in the direction of x
     */
    pub fn  next_after(&self,  x: &Dfp) -> Dfp  {
        // make sure we don't mix number with different precision
        if field.get_radix_digits() != x.field.get_radix_digits() {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INVALID);
             let result: Dfp = new_instance(get_zero());
            result.nans = QNAN;
            return dotrap(DfpField::FLAG_INVALID, NEXT_AFTER_TRAP, x, result);
        }
        // if this is greater than x
         let up: bool = false;
        if self.less_than(x) {
            up = true;
        }
        if compare(self, x) == 0 {
            return new_instance(x);
        }
        if less_than(get_zero()) {
            up = !up;
        }
         let inc: Dfp;
         let result: Dfp;
        if up {
            inc = new_instance(get_one());
            inc.exp = self.exp - mant.length + 1;
            inc.sign = self.sign;
            if self.equals(get_zero()) {
                inc.exp = MIN_EXP - mant.length;
            }
            result = add(inc);
        } else {
            inc = new_instance(get_one());
            inc.exp = self.exp;
            inc.sign = self.sign;
            if self.equals(inc) {
                inc.exp = self.exp - mant.length;
            } else {
                inc.exp = self.exp - mant.length + 1;
            }
            if self.equals(get_zero()) {
                inc.exp = MIN_EXP - mant.length;
            }
            result = self.subtract(inc);
        }
        if result.classify() == INFINITE && self.classify() != INFINITE {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
            result = dotrap(DfpField::FLAG_INEXACT, NEXT_AFTER_TRAP, x, result);
        }
        if result.equals(get_zero()) && self.equals(get_zero()) == false {
            field.set_i_e_e_e_flags_bits(DfpField::FLAG_INEXACT);
            result = dotrap(DfpField::FLAG_INEXACT, NEXT_AFTER_TRAP, x, result);
        }
        return result;
    }

    /** Convert the instance into a double.
     * @return a double approximating the instance
     * @see #toSplitDouble()
     */
    pub fn  to_double(&self) -> f64  {
        if is_infinite() {
            if less_than(get_zero()) {
                return Double::NEGATIVE_INFINITY;
            } else {
                return Double::POSITIVE_INFINITY;
            }
        }
        if is_na_n() {
            return Double::NaN;
        }
         let y: Dfp = self;
         let negate: bool = false;
         let cmp0: i32 = compare(self, get_zero());
        if cmp0 == 0 {
            return  if sign < 0 { -0.0 } else { 0.0 };
        } else if cmp0 < 0 {
            y = negate();
            negate = true;
        }
        /* Find the exponent, first estimate by integer log10, then adjust.
         Should be faster than doing a natural logarithm.  */
         let exponent: i32 = (y.int_log10() * 3.32) as i32;
        if exponent < 0 {
            exponent-= 1;
        }
         let temp_dfp: Dfp = DfpMath.pow(get_two(), exponent);
        while temp_dfp.less_than(y) || temp_dfp.equals(y) {
            temp_dfp = temp_dfp.multiply(2);
            exponent+= 1;
        }
        exponent-= 1;
        /* We have the exponent, now work on the mantissa */
        y = y.divide(DfpMath.pow(get_two(), exponent));
        if exponent > -1023 {
            y = y.subtract(get_one());
        }
        if exponent < -1074 {
            return 0;
        }
        if exponent > 1023 {
            return  if negate { Double::NEGATIVE_INFINITY } else { Double::POSITIVE_INFINITY };
        }
        y = y.multiply(new_instance(4503599627370496)).rint();
         let str: String = y.to_string();
        str = str.substring(0, str.length() - 1);
         let mantissa: i64 = Long.parse_long(str);
        if mantissa == 4503599627370496 {
            // Handle special case where we round up to next power of two
            mantissa = 0;
            exponent+= 1;
        }
        /* Its going to be subnormal, so make adjustments */
        if exponent <= -1023 {
            exponent-= 1;
        }
        while exponent < -1023 {
            exponent+= 1;
            mantissa >>= /* >>>= */ 1;
        }
         let bits: i64 = mantissa | ((exponent + 1023) << 52);
         let x: f64 = Double.long_bits_to_double(bits);
        if negate {
            x = -x;
        }
        return x;
    }

    /** Convert the instance into a split double.
     * @return an array of two doubles which sum represent the instance
     * @see #toDouble()
     */
    pub fn  to_split_double(&self) -> f64[]  {
         let split: [f64; 2] = [0.0; 2];
         let mask: i64 = 0xffffffffc0000000;
        split[0] = Double.long_bits_to_double(Double.double_to_long_bits(to_double()) & mask);
        split[1] = subtract(new_instance(split[0])).to_double();
        return split;
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  get_real(&self) -> f64  {
        return to_double();
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  add(&self,  a: f64) -> Dfp  {
        return add(new_instance(a));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  subtract(&self,  a: f64) -> Dfp  {
        return subtract(new_instance(a));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  multiply(&self,  a: f64) -> Dfp  {
        return multiply(new_instance(a));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  divide(&self,  a: f64) -> Dfp  {
        return divide(new_instance(a));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  remainder(&self,  a: f64) -> Dfp  {
        return remainder(new_instance(a));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  round(&self) -> i64  {
        return FastMath.round(to_double());
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  signum(&self) -> Dfp  {
        if is_na_n() || is_zero() {
            return self;
        } else {
            return new_instance( if sign > 0 { 1 } else { -1 });
        }
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  copy_sign(&self,  s: &Dfp) -> Dfp  {
        if (sign >= 0 && s.sign >= 0) || (sign < 0 && s.sign < 0) {
            // Sign is currently OK
            return self;
        }
        // flip sign
        return negate();
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  copy_sign(&self,  s: f64) -> Dfp  {
         let sb: i64 = Double.double_to_long_bits(s);
        if (sign >= 0 && sb >= 0) || (sign < 0 && sb < 0) {
            // Sign is currently OK
            return self;
        }
        // flip sign
        return negate();
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  scalb(&self,  n: i32) -> Dfp  {
        return multiply(DfpMath.pow(get_two(), n));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  hypot(&self,  y: &Dfp) -> Dfp  {
        return multiply(self).add(y.multiply(y)).sqrt();
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  cbrt(&self) -> Dfp  {
        return root_n(3);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  root_n(&self,  n: i32) -> Dfp  {
        return  if (sign >= 0) { DfpMath.pow(self, get_one().divide(n)) } else { DfpMath.pow(negate(), get_one().divide(n)).negate() };
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  pow(&self,  p: f64) -> Dfp  {
        return DfpMath.pow(self, new_instance(p));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  pow(&self,  n: i32) -> Dfp  {
        return DfpMath.pow(self, n);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  pow(&self,  e: &Dfp) -> Dfp  {
        return DfpMath.pow(self, e);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  exp(&self) -> Dfp  {
        return DfpMath.exp(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  expm1(&self) -> Dfp  {
        return DfpMath.exp(self).subtract(get_one());
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  log(&self) -> Dfp  {
        return DfpMath.log(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  log1p(&self) -> Dfp  {
        return DfpMath.log(self.add(get_one()));
    }

    //  TODO: deactivate this implementation (and return type) in 4.0
    /** Get the exponent of the greatest power of 10 that is less than or equal to abs(this).
     *  @return integer base 10 logarithm
     *  @deprecated as of 3.2, replaced by {@link #intLog10()}, in 4.0 the return type
     *  will be changed to Dfp
     */
    pub fn  log10(&self) -> i32  {
        return int_log10();
    }

    //    TODO: activate this implementation (and return type) in 4.0
    //    /** {@inheritDoc}
    //     * @since 3.2
    //     */
    //    public Dfp log10() {
    //        return DfpMath.log(this).divide(DfpMath.log(newInstance(10)));
    //    }
    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  cos(&self) -> Dfp  {
        return DfpMath.cos(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  sin(&self) -> Dfp  {
        return DfpMath.sin(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  tan(&self) -> Dfp  {
        return DfpMath.tan(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  acos(&self) -> Dfp  {
        return DfpMath.acos(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  asin(&self) -> Dfp  {
        return DfpMath.asin(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  atan(&self) -> Dfp  {
        return DfpMath.atan(self);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  atan2(&self,  x: &Dfp) -> /*  throws DimensionMismatchException */Result<Dfp>   {
        // compute r = sqrt(x^2+y^2)
         let r: Dfp = x.multiply(x).add(multiply(self)).sqrt();
        if x.sign >= 0 {
            // compute atan2(y, x) = 2 atan(y / (r + x))
            return get_two().multiply(divide(r.add(x)).atan());
        } else {
            // compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
             let tmp: Dfp = get_two().multiply(divide(r.subtract(x)).atan());
             let pm_pi: Dfp = new_instance( if (tmp.sign <= 0) { -FastMath::PI } else { FastMath::PI });
            return pm_pi.subtract(tmp);
        }
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  cosh(&self) -> Dfp  {
        return DfpMath.exp(self).add(DfpMath.exp(negate())).divide(2);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  sinh(&self) -> Dfp  {
        return DfpMath.exp(self).subtract(DfpMath.exp(negate())).divide(2);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  tanh(&self) -> Dfp  {
         let e_plus: Dfp = DfpMath.exp(self);
         let e_minus: Dfp = DfpMath.exp(negate());
        return e_plus.subtract(e_minus).divide(e_plus.add(e_minus));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  acosh(&self) -> Dfp  {
        return multiply(self).subtract(get_one()).sqrt().add(self).log();
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  asinh(&self) -> Dfp  {
        return multiply(self).add(get_one()).sqrt().add(self).log();
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  atanh(&self) -> Dfp  {
        return get_one().add(self).divide(get_one().subtract(self)).log().divide(2);
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a: &Dfp[],  b: &Dfp[]) -> /*  throws DimensionMismatchException */Result<Dfp>   {
        if a.length != b.length {
            throw new DimensionMismatchException(a.length, b.length);
        }
         let r: Dfp = get_zero();
        for ( let i: i32 = 0; i < a.length; i+= 1) {
            r = r.add(a[i].multiply(b[i]));
        }
        return r;
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a: &f64[],  b: &Dfp[]) -> /*  throws DimensionMismatchException */Result<Dfp>   {
        if a.length != b.length {
            throw new DimensionMismatchException(a.length, b.length);
        }
         let r: Dfp = get_zero();
        for ( let i: i32 = 0; i < a.length; i+= 1) {
            r = r.add(b[i].multiply(a[i]));
        }
        return r;
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a1: &Dfp,  b1: &Dfp,  a2: &Dfp,  b2: &Dfp) -> Dfp  {
        return a1.multiply(b1).add(a2.multiply(b2));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a1: f64,  b1: &Dfp,  a2: f64,  b2: &Dfp) -> Dfp  {
        return b1.multiply(a1).add(b2.multiply(a2));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a1: &Dfp,  b1: &Dfp,  a2: &Dfp,  b2: &Dfp,  a3: &Dfp,  b3: &Dfp) -> Dfp  {
        return a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a1: f64,  b1: &Dfp,  a2: f64,  b2: &Dfp,  a3: f64,  b3: &Dfp) -> Dfp  {
        return b1.multiply(a1).add(b2.multiply(a2)).add(b3.multiply(a3));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a1: &Dfp,  b1: &Dfp,  a2: &Dfp,  b2: &Dfp,  a3: &Dfp,  b3: &Dfp,  a4: &Dfp,  b4: &Dfp) -> Dfp  {
        return a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3)).add(a4.multiply(b4));
    }

    /** {@inheritDoc}
     * @since 3.2
     */
    pub fn  linear_combination(&self,  a1: f64,  b1: &Dfp,  a2: f64,  b2: &Dfp,  a3: f64,  b3: &Dfp,  a4: f64,  b4: &Dfp) -> Dfp  {
        return b1.multiply(a1).add(b2.multiply(a2)).add(b3.multiply(a3)).add(b4.multiply(a4));
    }
}


                
                