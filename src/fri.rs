use std::{str::FromStr, vec};

use num_bigint::{BigInt, BigUint};

use crate::{fft::FFT, prime_field::PrimeField};

pub struct FRI {
    field: PrimeField,
    fft: FFT,
}

impl FRI {
    pub fn new() -> FRI {
        let field = PrimeField::new(
            BigInt::from_str(
                "18446744069414584321",
            )
            .unwrap(),
        );
        Self {
            fft: FFT::new(field.clone()),
            field: field,
        }
    }

    pub fn construct(&self) {
        let root_of_unity = BigInt::from_str(
            "1753635133440165772",
        )
        .unwrap();
        let polynomial_to_test = vec![
            BigInt::from(34),
            BigInt::from(2),
            BigInt::from(10),
            BigInt::from(17),
            BigInt::from(9),
        ];

        // Domain for FRI evaluation and Folding (L0)
        let domain = self.fft.build_root_of_unity(&root_of_unity);
        self.fold(domain, polynomial_to_test, BigInt::from(10));
    }

    fn fold(&self, domain: Vec<BigInt>, polynomial: Vec<BigInt>, x0: BigInt) {
        let new_domain = self.fold_domain(&domain);

        if new_domain.len() <= 2 {
            println!("Domain length reached 2 exiting");
            return;
        }

        let new_domain_ref: Vec<&BigInt> = new_domain.iter().collect();

        //Split the polynomial
        let even_part: Vec<&BigInt> = polynomial.iter().step_by(2).collect();
        let odd_part: Vec<&BigInt> = polynomial.iter().skip(1).step_by(2).collect();

        //Evaluate on this new domain so that we get the x and y values of next round
        let even_part_eval = self.fft.compute_fft(&even_part, &new_domain_ref);
        let mut odd_part_eval = self.fft.compute_fft(&odd_part, &new_domain_ref);

        let p0 = even_part_eval;
        let p1: Vec<BigInt> = odd_part_eval
            .iter_mut()
            .map(|v| self.field.mul(&x0, v))
            .collect();

        let y = self.field.add_polynomial(&p0, &p1);
        self.fold(new_domain, y, x0);
    }

    fn fold_domain(&self, domain: &[BigInt]) -> Vec<BigInt> {
        //FIXME: Hard coding the shrinking polynomial for now
        // y = x^2
        let new_domain: Vec<BigInt> = domain.iter().map(|d| self.field.mul(d, d)).collect();
        new_domain[..domain.len() / 2].to_vec()
    }
}

#[cfg(test)]
mod tests {
    use crate::fri::FRI;

    #[test]
    pub fn test_fri() {
        let fri = FRI::new();
        fri.construct();
    }
}
