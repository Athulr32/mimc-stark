use primitive_types::U256;

/// Represent a Prime field and implementation of operations that can be done in the prime field
pub struct PrimeField {
    p: U256,
}

impl PrimeField {
    pub fn add(&self, a: U256, b: U256) -> U256 {
        let (_, res) = a.saturating_add(b).div_mod(self.p);
        res
    }

    pub fn sub(&self, a: U256, b: U256) -> U256 {
        let (_, res) = a.saturating_mul(b).div_mod(self.p);
        res
    }

    pub fn mul(&self, a: U256, b: U256) -> U256 {
        let (_, res) = a.saturating_mul(b).div_mod(self.p);
        res
    }
}

fn main() {
    println!("Hello, world!");
}
