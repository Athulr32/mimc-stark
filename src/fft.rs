use num_bigint::BigInt;
use num_traits::ops::euclid::Euclid;

use crate::prime_field::PrimeField;

pub struct FFT {
    p: BigInt,
}

impl FFT {
    pub fn new(field: PrimeField) -> FFT {
        Self { p: field.get_p() }
    }

    /// Build all the powers of root of unity (1, w, w^1, w^2...w^n-1, 1)
    /// Used as X coordinate for FFT
    pub fn build_root_of_unity(&self, root_of_unity: &BigInt) -> Vec<BigInt> {
        let mut result = vec![BigInt::from(1), root_of_unity.clone()];

        // We continue until we get w^N = 1
        while result.last().unwrap() != &BigInt::from(1) {
            result.push((result.last().unwrap() * root_of_unity) % &self.p);
        }

        result
    }

    /// Compute the DFT for the given polynomials at the roots
    fn compute_dft(&self, values: &[&BigInt], roots: &[&BigInt]) -> Vec<BigInt> {
        assert_eq!(values.len(), roots.len());
        let mut output = Vec::new();
        let len = roots.len();
        for i in 0..len {
            let mut last = BigInt::ZERO;
            for j in 0..len {
                last += values[j] * roots[(i * j) % len];
            }
            output.push(last % &self.p);
        }

        output
    }

    /// Evaluate the polynomial at roots of unity by computing fft
    /// we split the polynomial into odd and even indexes and recursively compute
    /// f(x) = feven(x^2) + x * fodd(x^2)
    pub fn compute_fft(&self, values: &[&BigInt], roots_of_unity: &[&BigInt]) -> Vec<BigInt> {
        if values.len() == 4 {
            return self.compute_dft(values, &roots_of_unity);
        }

        // Split into odd and even indices values
        let odd_values: Vec<&BigInt> = values.iter().skip(1).step_by(2).copied().collect();
        let even_values: Vec<&BigInt> = values.iter().step_by(2).copied().collect();

        // We need roots of unity at x^2,x^4...
        // because we compute for feven(x^2) not feven(x)
        let required_roots_of_unity: Vec<&BigInt> =
            roots_of_unity.iter().step_by(2).copied().collect();

        let l = self.compute_fft(&even_values, &required_roots_of_unity);
        let r = self.compute_fft(&odd_values, &required_roots_of_unity);

        let mut output = vec![BigInt::ZERO; values.len()];

        // Butterfly
        for (i, (x, y)) in l.iter().zip(r).enumerate() {
            // x = f_even(ω^{2i}) ← DFT of even part at ω^{2i}
            // y = f_odd(ω^{2i})  ← DFT of odd part at ω^{2i}

            //f_odd(w) = f(w^w) * w
            let f_odd = y * roots_of_unity[i];

            // f(w) = f_even(w^2) + w * f_odd(w^2)
            let fw = x + &f_odd;

            // for k = i
            output[i] = fw % &self.p;

            //for k = i + n/2
            output[i + l.len()] = (x - &f_odd).rem_euclid(&self.p);
        }

        output
    }

    // Padding the values based on the length of roots it should have same size
    // the +1 after .len() is to account the last element of roots that is 1 which is not needed for computation
    pub fn pad(&self, roots: &[BigInt], values: &mut Vec<BigInt>) {
        if roots.len() > values.len() + 1 {
            let n = roots.len() - values.len() - 1;
            values.extend_from_slice(&vec![BigInt::ZERO; n]);
        }
    }

    /// Fast Fourier Transform
    pub fn fft(&self, mut values: Vec<BigInt>, root_of_unity: &BigInt) -> Vec<BigInt> {
        let roots = self.build_root_of_unity(root_of_unity);

        self.pad(&roots, &mut values);
        assert_eq!(roots.len(), values.len() + 1);

        let values: Vec<&BigInt> = values.iter().collect();
        let roots: Vec<&BigInt> = roots.iter().collect();

        return self.compute_fft(&values, &roots[..roots.len() - 1]);
    }

    /// Inverse Fast Fourier Transform
    pub fn inv_fft(&self, mut values: Vec<BigInt>, root_of_unity: &BigInt) -> Vec<BigInt> {
        let mut roots = self.build_root_of_unity(root_of_unity);

        self.pad(&roots, &mut values);

        let values: Vec<&BigInt> = values.iter().collect();
        roots.reverse();
        let roots: Vec<&BigInt> = roots[..roots.len() - 1].iter().collect();

        let n_inv = BigInt::from(values.len()).modpow(&(&self.p - 2), &self.p); //Modular invesrse of values.len()
        let mut inv_fft = self.compute_fft(&values, &roots);
        inv_fft.iter_mut().for_each(|inv| {
            *inv = (&*inv * &n_inv) % &self.p;
        });
        return inv_fft;
    }

    /// Multiply two polynomials
    pub fn multiply_polynomials(
        &self,
        mut a: Vec<BigInt>,
        mut b: Vec<BigInt>,
        root_of_unity: &BigInt,
    ) -> Vec<BigInt> {
        let roots = self.build_root_of_unity(root_of_unity);

        self.pad(&roots, &mut a);
        self.pad(&roots, &mut b);
        let roots_ref: Vec<&BigInt> = roots.iter().collect();
        let values_a: Vec<&BigInt> = a.iter().collect();
        let values_b: Vec<&BigInt> = b.iter().collect();

        let x1 = self.compute_fft(&values_a, &roots_ref[..roots_ref.len() - 1]);
        let x2 = self.compute_fft(&values_b, &roots_ref[..roots_ref.len() - 1]);

        let y: Vec<BigInt> = x1
            .iter()
            .zip(x2)
            .map(|(x1, x2)| (x1 * x2) % &self.p)
            .collect();

        self.inv_fft(y, root_of_unity)
    }
}

#[cfg(test)]
mod test {
    use num_bigint::{BigInt, BigUint};

    use crate::fft::FFT;

    #[test]
    fn test_root_of_unity() {
        let root_of_unity = BigInt::from(4);
        let p = BigInt::from(17);
        let fft = FFT { p };
        let result = fft.build_root_of_unity(&root_of_unity);
        assert_eq!(
            result,
            vec![
                BigInt::from(1),
                BigInt::from(4),
                BigInt::from(16),
                BigInt::from(13),
                BigInt::from(1)
            ]
        );
    }

    #[test]
    fn test_values_pad() {
        let root_of_unity = BigInt::from(4);
        let p = BigInt::from(17);
        let fft = FFT { p };
        let roots = fft.build_root_of_unity(&root_of_unity);
        let mut values = vec![BigInt::from(2), BigInt::from(43)];
        fft.pad(&roots, &mut values);
        assert_eq!(roots.len(), values.len() + 1);
    }

    #[test]
    fn test_fft() {
        let root_of_unity = BigInt::from(9);
        let poly = vec![
            BigInt::from(1),
            BigInt::from(2),
            BigInt::from(3),
            BigInt::from(4),
            BigInt::from(5),
            BigInt::from(6),
        ];

        let p = BigInt::from(17);
        let fft = FFT { p };
        let result = fft.fft(poly, &root_of_unity);
        assert_eq!(
            result,
            vec![
                BigInt::from(4),
                BigInt::from(8),
                BigInt::from(4),
                BigInt::from(1),
                BigInt::from(14),
                BigInt::from(11),
                BigInt::from(2),
                BigInt::from(15)
            ]
        );
    }

    #[test]
    fn test_mul_poly() {
        let a: Vec<BigInt> = vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)];

        let b: Vec<BigInt> = vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)];
        let root_of_unity = BigInt::from(9);
        let p = BigInt::from(17);

        let fft = FFT { p };
        let result = fft.multiply_polynomials(a, b, &root_of_unity);
        println!("{:?}", result);
    }
}
