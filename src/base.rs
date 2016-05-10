
pub trait FieldElement<T> {
	
	 /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    fn  add(&self,  a: f64) -> T ;

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    fn  subtract(&self,  a: f64) -> T ;

	fn negate(&self) -> T;
	
	fn multiply_by_integer(&self, n: i64) -> T;
	
	fn multiply(&self, a: T) -> T;
	
	fn divide(&self, a: T) -> Result<T, &'static str>;  // MathArithmeticsException
	  
	fn reciprocal(&self) -> Result<T, &'static str>;  // MathArithmeticsException
	
	// fn getField() -> &'static Field<T>;
	
}

/**
 * Interface representing a <a href="http://mathworld.wolfram.com/RealNumber.html">real</a>
 * <a href="http://mathworld.wolfram.com/Field.html">field</a>.
 * @param <T> the type of the field elements
 * @see FieldElement
 * @since 3.2
 */
pub trait RealFieldElement<T> : FieldElement<T> {

    /** Get the real value of the number.
     * @return real value
     */
    fn  get_real(&self) -> f64 ;

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    fn  add(&self,  a: f64) -> T ;

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    fn  subtract(&self,  a: f64) -> T ;

    /** '×' operator.
     * @param a right hand side parameter of the operator
     * @return this×a
     */
    fn  multiply(&self,  a: f64) -> T ;

    /** '÷' operator.
     * @param a right hand side parameter of the operator
     * @return this÷a
     */
    fn  divide(&self,  a: f64) -> T ;

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n × a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    fn  remainder_f64(&self,  a: f64) -> T ;

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n × a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     * @exception DimensionMismatchException if number of free parameters or orders are inconsistent
     */
    fn  remainder(&self,  a: &T) -> /*  throws DimensionMismatchException */Result<T, &'static str>  ;

    /** absolute value.
     * @return abs(this)
     */
    fn  abs(&self) -> T ;

    /** Get the smallest whole number larger than instance.
     * @return ceil(this)
     */
    fn  ceil(&self) -> T ;

    /** Get the largest whole number smaller than instance.
     * @return floor(this)
     */
    fn  floor(&self) -> T ;

    /** Get the whole number that is the nearest to the instance, or the even one if x is exactly half way between two integers.
     * @return a double number r such that r is an integer r - 0.5 ≤ this ≤ r + 0.5
     */
    fn  rint(&self) -> T ;

    /** Get the closest long to instance value.
     * @return closest long to {@link #getReal()}
     */
    fn  round(&self) -> i64 ;

    /** Compute the signum of the instance.
     * The signum is -1 for negative numbers, +1 for positive numbers and 0 otherwise
     * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
     */
    fn  signum(&self) -> T ;

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    fn  copy_sign(&self,  sign: &T) -> T ;

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    fn  copy_sign_f64(&self,  sign: f64) -> T ;

    /**
     * Multiply the instance by a power of 2.
     * @param n power of 2
     * @return this × 2<sup>n</sup>
     */
    fn  scalb(&self,  n: i32) -> T ;

    /**
     * Returns the hypotenuse of a triangle with sides {@code this} and {@code y}
     * - sqrt(<i>this</i><sup>2</sup> +<i>y</i><sup>2</sup>)
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * @param y a value
     * @return sqrt(<i>this</i><sup>2</sup> +<i>y</i><sup>2</sup>)
     * @exception DimensionMismatchException if number of free parameters or orders are inconsistent
     */
    fn  hypot(&self,  y: &T) -> /*  throws DimensionMismatchException */Result<T, &'static str>  ;

    /** {@inheritDoc} */
    fn  reciprocal(&self) -> T ;

    /** Square root.
     * @return square root of the instance
     */
    fn  sqrt(&self) -> T ;

    /** Cubic root.
     * @return cubic root of the instance
     */
    fn  cbrt(&self) -> T ;

    /** N<sup>th</sup> root.
     * @param n order of the root
     * @return n<sup>th</sup> root of the instance
     */
    fn  root_n(&self,  n: i32) -> T ;

    /** Power operation.
     * @param p power to apply
     * @return this<sup>p</sup>
     */
    fn  pow_f64(&self,  p: f64) -> T ;

    /** Integer power operation.
     * @param n power to apply
     * @return this<sup>n</sup>
     */
    fn  pow_i32(&self,  n: i32) -> T ;

    /** Power operation.
     * @param e exponent
     * @return this<sup>e</sup>
     * @exception DimensionMismatchException if number of free parameters or orders are inconsistent
     */
    fn  pow(&self,  e: &T) -> /*  throws DimensionMismatchException */Result<T, &'static str>  ;

    /** Exponential.
     * @return exponential of the instance
     */
    fn  exp(&self) -> T ;

    /** Exponential minus 1.
     * @return exponential minus one of the instance
     */
    fn  expm1(&self) -> T ;

    /** Natural logarithm.
     * @return logarithm of the instance
     */
    fn  log(&self) -> T ;

    /** Shifted natural logarithm.
     * @return logarithm of one plus the instance
     */
    fn  log1p(&self) -> T ;

    //    TODO: add this method in 4.0, as it is not possible to do it in 3.2
    //          due to incompatibility of the return type in the Dfp class
    //    /** Base 10 logarithm.
    //     * @return base 10 logarithm of the instance
    //     */
    //    T log10();
    /** Cosine operation.
     * @return cos(this)
     */
    fn  cos(&self) -> T ;

    /** Sine operation.
     * @return sin(this)
     */
    fn  sin(&self) -> T ;

    /** Tangent operation.
     * @return tan(this)
     */
    fn  tan(&self) -> T ;

    /** Arc cosine operation.
     * @return acos(this)
     */
    fn  acos(&self) -> T ;

    /** Arc sine operation.
     * @return asin(this)
     */
    fn  asin(&self) -> T ;

    /** Arc tangent operation.
     * @return atan(this)
     */
    fn  atan(&self) -> T ;

    /** Two arguments arc tangent operation.
     * @param x second argument of the arc tangent
     * @return atan2(this, x)
     * @exception DimensionMismatchException if number of free parameters or orders are inconsistent
     */
    fn  atan2(&self,  x: &T) -> /*  throws DimensionMismatchException */Result<T, &'static str>  ;

    /** Hyperbolic cosine operation.
     * @return cosh(this)
     */
    fn  cosh(&self) -> T ;

    /** Hyperbolic sine operation.
     * @return sinh(this)
     */
    fn  sinh(&self) -> T ;

    /** Hyperbolic tangent operation.
     * @return tanh(this)
     */
    fn  tanh(&self) -> T ;

    /** Inverse hyperbolic cosine operation.
     * @return acosh(this)
     */
    fn  acosh(&self) -> T ;

    /** Inverse hyperbolic sine operation.
     * @return asin(this)
     */
    fn  asinh(&self) -> T ;

    /** Inverse hyperbolic  tangent operation.
     * @return atanh(this)
     */
    fn  atanh(&self) -> T ;

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>Σ<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws DimensionMismatchException if arrays dimensions don't match
     * @since 3.2
     */
    fn  linear_combination_vec(&self,  a: &Vec<T>,  b: &Vec<T>) -> /*  throws DimensionMismatchException */Result<T, &'static str>  ;

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>Σ<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws DimensionMismatchException if arrays dimensions don't match
     * @since 3.2
     */
    fn  linear_combination_f64_vec(&self,  a: &Vec<f64>,  b: &Vec<T>) -> /*  throws DimensionMismatchException */Result<T, &'static str>  ;

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @return a<sub>1</sub>×b<sub>1</sub> +
     * a<sub>2</sub>×b<sub>2</sub>
     * @see #linearCombination(Object, Object, Object, Object, Object, Object)
     * @see #linearCombination(Object, Object, Object, Object, Object, Object, Object, Object)
     * @since 3.2
     */
    fn  linear_combination_2_factors(&self,  a1: &T,  b1: &T,  a2: &T,  b2: &T) -> T ;

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @return a<sub>1</sub>×b<sub>1</sub> +
     * a<sub>2</sub>×b<sub>2</sub>
     * @see #linearCombination(double, Object, double, Object, double, Object)
     * @see #linearCombination(double, Object, double, Object, double, Object, double, Object)
     * @since 3.2
     */
    fn  linear_combination_f64_2_factors(&self,  a1: f64,  b1: &T,  a2: f64,  b2: &T) -> T ;

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>×b<sub>1</sub> +
     * a<sub>2</sub>×b<sub>2</sub> + a<sub>3</sub>×b<sub>3</sub>
     * @see #linearCombination(Object, Object, Object, Object)
     * @see #linearCombination(Object, Object, Object, Object, Object, Object, Object, Object)
     * @since 3.2
     */
    fn  linear_combination_3_factors(&self,  a1: &T,  b1: &T,  a2: &T,  b2: &T,  a3: &T,  b3: &T) -> T ;

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>×b<sub>1</sub> +
     * a<sub>2</sub>×b<sub>2</sub> + a<sub>3</sub>×b<sub>3</sub>
     * @see #linearCombination(double, Object, double, Object)
     * @see #linearCombination(double, Object, double, Object, double, Object, double, Object)
     * @since 3.2
     */
    fn  linear_combination_f64_3_factors(&self,  a1: f64,  b1: &T,  a2: f64,  b2: &T,  a3: f64,  b3: &T) -> T ;

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @param a4 first factor of the third term
     * @param b4 second factor of the third term
     * @return a<sub>1</sub>×b<sub>1</sub> +
     * a<sub>2</sub>×b<sub>2</sub> + a<sub>3</sub>×b<sub>3</sub> +
     * a<sub>4</sub>×b<sub>4</sub>
     * @see #linearCombination(Object, Object, Object, Object)
     * @see #linearCombination(Object, Object, Object, Object, Object, Object)
     * @since 3.2
     */
    fn  linear_combination_4_factors(&self,  a1: &T,  b1: &T,  a2: &T,  b2: &T,  a3: &T,  b3: &T,  a4: &T,  b4: &T) -> T ;

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @param a4 first factor of the third term
     * @param b4 second factor of the third term
     * @return a<sub>1</sub>×b<sub>1</sub> +
     * a<sub>2</sub>×b<sub>2</sub> + a<sub>3</sub>×b<sub>3</sub> +
     * a<sub>4</sub>×b<sub>4</sub>
     * @see #linearCombination(double, Object, double, Object)
     * @see #linearCombination(double, Object, double, Object, double, Object)
     * @since 3.2
     */
    fn  linear_combination_f64_4_factors(&self,  a1: f64,  b1: &T,  a2: f64,  b2: &T,  a3: f64,  b3: &T,  a4: f64,  b4: &T) -> T ;
}



pub trait Field<T> {
	
	fn get_zero() -> T;
	
	fn get_one() -> T;
	
}