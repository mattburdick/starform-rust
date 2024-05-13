// src/types.rs

use std::fmt;

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum MassType {
    Star,
    Planet,
    Moon,
    GasGiant,
}

impl fmt::Display for MassType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            MassType::Star => write!(f, "Star"),
            MassType::Planet => write!(f, "Planet"),
            MassType::Moon => write!(f, "Moon"),
            MassType::GasGiant => write!(f, "Gas Giant"),
        }
    }
}
