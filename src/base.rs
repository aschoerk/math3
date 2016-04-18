
pub trait FieldElement<T> {
	
	fn add(a: T) -> T;
	
	fn subtract(a: T) ->T;
	
	fn negate() -> T;
	
	fn multiply_by_integer(n: i64) -> T;
	
	fn multiply(a: T) -> T;
	
	fn divide(a: T) -> Result<T, &'static str>;  // MathArithmeticsException
	  
	fn reciprocal() -> Result<T, &'static str>;  // MathArithmeticsException
	
	// fn getField() -> &'static Field<T>;
	
}


pub trait Field<T> {
	fn get_zero() -> T;
	
	fn get_one() -> T;
	
}