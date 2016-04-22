#![allow(dead_code)] 
use std::f64;
use rsutils::*;

/** Archimede's constant PI, ratio of circle circumference to diameter. */
pub const PI: f64 = 105414357.0 / 33554432.0 + 1.984187159361080883e-9;
/** Napier's constant e, base of the natural logarithm. */
pub const E: f64 = 2850325.0 / 1048576.0 + 8.254840070411028747e-8;

/** Index of exp(0) in the array of integer exponentials. */
const EXP_INT_TABLE_MAX_INDEX: i16 = 750;
/** Length of the array of integer exponentials. */
const EXP_INT_TABLE_LEN: i16 = EXP_INT_TABLE_MAX_INDEX * 2;
/** Logarithm table length. */
const LN_MANT_LEN: i16 = 1024;
/** Exponential fractions table length. */
const EXP_FRAC_TABLE_LEN: i16 = 1025; // 0, 1/1024, ... 1024/1024

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
    const SINE_TABLE_LEN: i16 = 14;

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
    const  RECIP_2PI: [i64; 18] = [
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
    const PI_O_4_BITS: [i64; 2] = [
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
    const HEX_40000000: i64 = 0x40000000; // 1073741824L

    /** Mask used to clear low order 30 bits */
    const MASK_30BITS: i64 = -1 - (HEX_40000000 -1); // 0xFFFFFFFFC0000000L;

    /** Mask used to clear the non-sign part of an int. */
    const MASK_NON_SIGN_INT: i32 = 0x7fffffff;

    /** Mask used to clear the non-sign part of a long. */
    const MASK_NON_SIGN_LONG: i64 = 0x7fffffffffffffff;

    /** Mask used to extract exponent from double bits. */
    const MASK_DOUBLE_EXPONENT: i64 = 0x7ff0000000000000;

    /** Mask used to extract mantissa from double bits. */
    const MASK_DOUBLE_MANTISSA: i64 = 0x000fffffffffffff;

    /** Mask used to add implicit high order bit for normalized double. */
    const IMPLICIT_HIGH_BIT: i64 = 0x0010000000000000;

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

 
pub fn abs(x: f64) -> f64 {
	long_bits_to_double(&(MASK_NON_SIGN_LONG & double_to_raw_long_bits(&x)))
}


struct FastMath {
	
}