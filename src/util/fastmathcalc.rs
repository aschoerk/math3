#![allow(dead_code)] 



        /**
     * 0x40000000 - used to split a double into two parts, both with the low order bits cleared.
     * Equivalent to 2^30.
     */
    // 1073741824L
     const HEX_40000000: i64 = 0x40000000;

    /** Factorial table, for Taylor series expansions. 0!, 1!, 2!, ... 19! */
     const FACT : [f64; 20] = [// 0
    1.0, // 1
    1.0, // 2
    2.0, // 3
    6.0, // 4
    24.0, // 5
    120.0, // 6
    720.0, // 7
    5040.0, // 8
    40320.0, // 9
    362880.0, // 10
    3628800.0, // 11
    39916800.0, // 12
    479001600.0, // 13
    6227020800.0, // 14
    87178291200.0, // 15
    1307674368000.0, // 16
    20922789888000.0, // 17
    355687428096000.0, // 18
    6402373705728000.0, // 19
    121645100408832000.0, ]
    ;

    /** Coefficients for slowLog. */
     const LN_SPLIT_COEF: [[f64; 2]; 16] = [[2.0, 0.0, ]
    , [0.6666666269302368, 3.9736429850260626E-8, ]
    , [0.3999999761581421, 2.3841857910019882E-8, ]
    , [0.2857142686843872, 1.7029898543501842E-8, ]
    , [0.2222222089767456, 1.3245471311735498E-8, ]
    , [0.1818181574344635, 2.4384203044354907E-8, ]
    , [0.1538461446762085, 9.140260083262505E-9, ]
    , [0.13333332538604736, 9.220590270857665E-9, ]
    , [0.11764700710773468, 1.2393345855018391E-8, ]
    , [0.10526403784751892, 8.251545029714408E-9, ]
    , [0.0952233225107193, 1.2675934823758863E-8, ]
    , [0.08713622391223907, 1.1430250008909141E-8, ]
    , [0.07842259109020233, 2.404307984052299E-9, ]
    , [0.08371849358081818, 1.176342548272881E-8, ]
    , [0.030589580535888672, 1.2958646899018938E-9, ]
    , [0.14982303977012634, 1.225743062930824E-8, ]
    , ]
    ;

    /** Table start declaration. */
     const TABLE_START_DECL: &'static str = "    {";

    /** Table end declaration. */
     const TABLE_END_DECL: &'static str = "    };";
     
     struct FastMathCalc {
     }
     
     /**
     *  For x between 0 and 1, returns exp(x), uses extended precision
     *  @param x argument of exponential
     *  @param result placeholder where to place exp(x) split in two terms
     *  for extra precision (i.e. exp(x) = result[0] + result[1]
     *  @return exp(x)
     */
    pub fn  slowexp( x: f64) -> (f64, [f64; 2])  {
    	let mut result: [f64; 2] = [0.0; 2];
         let mut xs: [f64; 2];
         let mut ys: [f64; 2] = [0.0; 2];
         let mut facts: [f64; 2];
        xs = split(x);
        ys[1] = 0.0;
        ys[0] = 0.0;
        for i in (0..FACT.len()).rev() {
            let mut tmp = split_mult(xs, ys);
            ys[0] = tmp[0];
            ys[1] = tmp[1];
            tmp = split(FACT[i]);
            facts = split_reciprocal(tmp);
            tmp = split_add(ys, facts);
            ys[0] = tmp[0];
            ys[1] = tmp[1];
        }
        result[0] = ys[0];
        result[1] = ys[1];
        
        return (ys[0] + ys[1], result);
    }
     
      /** Compute split[0], split[1] such that their sum is equal to d,
     * and split[0] has its 30 least significant bits as zero.
     * @param d number to split
     * @param split placeholder where to place the result
     */
    fn  split( d: f64) -> [f64;2]  {
    	
        if d < 8e298 && d > -8e298 {
            let a: f64 = d * HEX_40000000 as f64;
            [(d + a) - a, d - ((d + a) - a)]
        } else {
             let a: f64 = d * 9.31322574615478515625E-10;
            let tmp = (d + a - d) * HEX_40000000 as f64;
            [tmp, d - tmp]
        }
    }

    /** Recompute a split.
     * @param a input/out array containing the split, changed
     * on output
     */
    fn  resplit( a: [f64;2]) -> [f64;2]  {
         let c: f64 = a[0] + a[1];
         let d: f64 = -(c - a[0] - a[1]);
        if c < 8e298 && c > -8e298 {
            // MAGIC NUMBER
            let z: f64 = c * HEX_40000000 as f64;
             
            let tmp = (c + z) - z;
            [tmp, c - tmp + d]
        } else {
             let z: f64 = c * 9.31322574615478515625E-10;
            let tmp = (c + z - c) * HEX_40000000 as f64;
            [tmp, c - tmp + d]
        }
    }

    /** Multiply two numbers in split form.
     * @param a first term of multiplication
     * @param b second term of multiplication
     * @param ans placeholder where to put the result
     */
    fn  split_mult( a: [f64; 2],  b: [f64; 2]) -> [f64; 2]  {
        let ans0 = a[0] * b[0];
        let ans1 = a[0] * b[1] + a[1] * b[0] + a[1] * b[1];
        /* Resplit */
        resplit([ans0, ans1])
    }

    /** Add two numbers in split form.
     * @param a first term of addition
     * @param b second term of addition
     * @param ans placeholder where to put the result
     */
    fn  split_add( a: [f64; 2],  b: [f64; 2]) -> [f64; 2]  {
        let ans0 = a[0] + b[0];
        let ans1 = a[1] + b[1];
        resplit([ans0, ans1])
    }

    /** Compute the reciprocal of in.  Use the following algorithm.
     *  in = c + d.
     *  want to find x + y such that x+y = 1/(c+d) and x is much
     *  larger than y and x has several zero bits on the right.
     *
     *  Set b = 1/(2^22),  a = 1 - b.  Thus (a+b) = 1.
     *  Use following identity to compute (a+b)/(c+d)
     *
     *  (a+b)/(c+d)  =   a/c   +    (bc - ad) / (c^2 + cd)
     *  set x = a/c  and y = (bc - ad) / (c^2 + cd)
     *  This will be close to the right answer, but there will be
     *  some rounding in the calculation of X.  So by carefully
     *  computing 1 - (c+d)(x+y) we can compute an error and
     *  add that back in.   This is done carefully so that terms
     *  of similar size are subtracted first.
     *  @param in initial number, in split form
     *  @param result placeholder where to put the result
     */
    pub fn  split_reciprocal( inp_p: [f64; 2]) -> [f64; 2]  {
    	let mut inp = inp_p;
         let b: f64 = 1.0 / 4194304.0;
         let a: f64 = 1.0 - b;
        if inp[0] == 0.0 {
            inp[0] = inp[1];
            inp[1] = 0.0;
        }
        let result0 = a / inp[0];
        let mut result1 = (b * inp[0] - a * inp[1]) / (inp[0] * inp[0] + inp[0] * inp[1]);
        if result1 != result1 {
            // can happen if result[1] is NAN
            result1 = 0.0;
        }
        /* Resplit */
        let mut result = resplit([result0, result1]);
        for i in 0..2 {
            /* this may be overkill, probably once is enough */
            let mut err: f64 = 1.0 - result[0] * inp[0] - result[0] * inp[1] - result[1] * inp[0] - result[1] * inp[1];
            /*err = 1.0 - err; */
            err *= result[0] + result[1];
            /*printf("err = %16e\n", err); */
            result[1] += err;
        }
        result
    }
    
     /** Compute (a[0] + a[1]) * (b[0] + b[1]) in extended precision.
     * @param a first term of the multiplication
     * @param b second term of the multiplication
     * @param result placeholder where to put the result
     */
    fn  quad_mult( a: [f64; 2],  b: [f64; 2]) -> [f64; 2]  {
    	let mut result: [f64; 2] = [0.0; 2];
         let mut xs: [f64; 2];
         let mut ys: [f64; 2];
         let mut zs: [f64; 2];
        /* a[0] * b[0] */
        xs = split(a[0]);
        ys = split(b[0]);
        zs = split_mult(xs, ys);
        result[0] = zs[0];
        result[1] = zs[1];
        /* a[0] * b[1] */
        ys = split(b[1]);
        zs = split_mult(xs, ys);
         let mut tmp: f64 = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;
        /* a[1] * b[0] */
        xs = split(a[1]);
        ys = split(b[0]);
        zs = split_mult(xs, ys);
        tmp = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;
        /* a[1] * b[0] */
        xs = split(a[1]);
        ys = split(b[1]);
        zs = split_mult(xs, ys);
        tmp = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;
        result
    }

    
    /** Compute exp(p) for a integer p in extended precision.
     * @param p integer whose exponential is requested
     * @param result placeholder where to put the result in extended precision
     * @return exp(p) in standard precision (equal to result[0] + result[1])
     */
    pub fn  expint( p_p: i32) -> (f64, [f64; 2])  {
    	let mut p = p_p;
    	let mut result: [f64;2] = [0.0; 2];
        //double x = M_E;
         let mut xs: [f64; 2] = [0.0; 2];
         let mut bs: [f64; 2] = [0.0; 2];
         let mut ys: [f64; 2] = [0.0; 2];
        //split(x, xs);
        //xs[1] = (double)(2.7182818284590452353602874713526625L - xs[0]);
        //xs[0] = 2.71827697753906250000;
        //xs[1] = 4.85091998273542816811e-06;
        //xs[0] = Double.longBitsToDouble(0x4005bf0800000000L);
        //xs[1] = Double.longBitsToDouble(0x3ed458a2bb4a9b00L);
        /* E */
        xs[0] = 2.718281828459045;
        xs[1] = 1.4456468917292502E-16;
        ys = split(1.0);
        while p > 0 {
            if (p & 1) != 0 {
                let tmp = quad_mult(ys, xs);
                ys[0] = tmp[0];
                ys[1] = tmp[1];
            }
            let tmp = quad_mult(xs, xs);
            xs[0] = tmp[0];
            xs[1] = tmp[1];
            p >>= 1;
        }
        result[0] = ys[0];
        result[1] = ys[1];
        result = resplit(result);
        return (ys[0] + ys[1], result);
    }

     
     /** xi in the range of [1, 2].
     *                                3        5        7
     *      x+1           /          x        x        x          \
     *  ln ----- =   2 *  |  x  +   ----  +  ----  +  ---- + ...  |
     *      1-x           \          3        5        7          /
     *
     * So, compute a Remez approximation of the following function
     *
     *  ln ((sqrt(x)+1)/(1-sqrt(x)))  /  x
     *
     * This will be an even function with only positive coefficents.
     * x is in the range [0 - 1/3].
     *
     * Transform xi for input to the above function by setting
     * x = (xi-1)/(xi+1).   Input to the polynomial is x^2, then
     * the result is multiplied by x.
     * @param xi number from which log is requested
     * @return log(xi)
     */
    pub fn  slow_log( xi: f64) -> [f64; 2]  {
        let mut x = split(xi);
        /* Set X = (x-1)/(x+1) */
        x[0] += 1.0;
        x = resplit(x);
        let mut a = split_reciprocal(x);
        x[0] -= 2.0;
        x = resplit(x);
        let mut y = split_mult(x, a);
        x[0] = y[0];
        x[1] = y[1];
        /* Square X -> X2*/
        let x2 = split_mult(x, x);
        //x[0] -= 1.0;
        //resplit(x);
        y[0] = LN_SPLIT_COEF[16 - 1][0];
        y[1] = LN_SPLIT_COEF[16 - 1][1];
        for i in (0..14).rev() {
            a = split_mult(y, x2);
            y[0] = a[0];
            y[1] = a[1];
            a = split_add(y, LN_SPLIT_COEF[i]);
            y[0] = a[0];
            y[1] = a[1];
        }
        a = split_mult(y, x);
        y[0] = a[0];
        y[1] = a[1];
        return y;
    }
     



