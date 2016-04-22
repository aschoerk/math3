#![allow(dead_code)] 
#[cfg(test)] extern crate hamcrest;

#[cfg(test)] #[macro_use] pub mod assert;

pub mod base;
pub mod fastmath;
// pub mod complex;
pub mod rsutils;
pub mod precision;


static a:i64 = (0x28be60db << 32) | 0x9391054a;

static b: f64 = 1.0 / 3.0;