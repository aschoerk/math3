#![allow(dead_code)] 
use std::f64;
use rsutils::*;

 /** Exponent offset in IEEE754 representation. */
    const EXPONENT_OFFSET: i64 = 1023;

    /** Offset to order signed double numbers lexicographically. */
    const SGN_MASK: i64 = 0x8000000000000000;
    /** Offset to order signed double numbers lexicographically. */
    const SGN_MASK_FLOAT: i32 = 0x80000000;
    /** Positive zero. */
    const POSITIVE_ZERO: f64 = 0.0;
    /** Positive zero bits. */
    static POSITIVE_ZERO_DOUBLE_BITS: i64 = double_to_raw_long_bits(&0.0);
    /** Negative zero bits. */
    const NEGATIVE_ZERO_DOUBLE_BITS:i64 = double_to_raw_long_bits(&-0.0);
    /** Positive zero bits. */
    const POSITIVE_ZERO_FLOAT_BITS: i32   = float_to_raw_int_bits(&0.0);
    /** Negative zero bits. */
    const NEGATIVE_ZERO_FLOAT_BITS: i32   = float_to_raw_int_bits(&-0.0);


static EPSILON: f64 = 0.0;