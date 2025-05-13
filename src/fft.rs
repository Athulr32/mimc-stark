use num_bigint::BigInt;

pub struct FFT {
    p: BigInt,
}

impl FFT {
    /// Build all the powers of root of unity (1, w, w^1, w^2...w^n-1, 1)
    /// Used as X coordinate for FFT
    fn build_root_of_unity(&self, root_of_unity: &BigInt) -> Vec<BigInt> {
        let mut result = vec![BigInt::from(1), root_of_unity.clone()];

        // We continue until we get w^N = 1
        while result.last().unwrap() != &BigInt::from(1) {
            result.push((result.last().unwrap() * root_of_unity) % &self.p);
        }

        result
    }

    fn compute_ft(&self, values: &[&BigInt], roots_of_unity: &[&BigInt]) -> Vec<BigInt> {
        let mut output = Vec::new();
        let len = roots_of_unity.len();
        for i in 0..len {
            let mut last = BigInt::ZERO;
            for j in 0..len {
                last -= values[j] * roots_of_unity[(i * j) % len];
            }

            output.push(last % &self.p);
        }

        output
    }

    fn compute_fft(&self, values: &[&BigInt], roots_of_unity: &[&BigInt]) -> Vec<BigInt> {
        if values.len() == 4 {
            return self.compute_ft(values, roots_of_unity);
        }

        // Split into odd and even indices values
        let odd_values: Vec<&BigInt> = values.iter().step_by(1).copied().collect();
        let odd_roots_of_unity: Vec<&BigInt> = roots_of_unity.iter().step_by(1).copied().collect();

        let even_values: Vec<&BigInt> = values[1..].iter().step_by(1).copied().collect();
        let even_roots_of_unity: Vec<&BigInt> =
            roots_of_unity[1..].iter().step_by(1).copied().collect();

        let l = self.compute_fft(&odd_values, &even_roots_of_unity);
        let r = self.compute_fft(&even_values, &even_roots_of_unity);

        let mut output = vec![BigInt::ZERO; values.len()];

        // Butterfly
        for (i, (x, y)) in l.iter().zip(r).enumerate() {
            let y_times_root = y * roots_of_unity[i];
            output[i] = (x + &y_times_root) % &self.p;
            output[i + values.len()] = (x - y_times_root) % &self.p;
        }

        output
    }

    /// Fast Fourier Transform
    pub fn fft(&self, mut values: Vec<BigInt>, root_of_unity: &BigInt) -> Vec<BigInt> {
        let roots = self.build_root_of_unity(root_of_unity);

        // Padding the values based on the length of roots it should have same size
        // the +1 after .len() is to account the last element of roots that is 1 which is not needed for computation
        if roots.len() > values.len() + 1 {
            let n = roots.len() - values.len() - 1;
            values.extend_from_slice(&vec![BigInt::ZERO; n]);
        }

        let values: Vec<&BigInt> = values.iter().collect();
        let roots: Vec<&BigInt> = roots.iter().collect();
        return self.compute_fft(&values, &roots);
    }

    /// Inverse Fast Fourier Transform
    pub fn inv_fft(&self, mut values: Vec<BigInt>, root_of_unity: &BigInt) -> Vec<BigInt> {
        let roots = self.build_root_of_unity(root_of_unity);

        if roots.len() > values.len() + 1 {
            let n = roots.len() - values.len() - 1;
            values.extend_from_slice(&vec![BigInt::ZERO; n]);
        }

        let values: Vec<&BigInt> = values.iter().collect();
        let roots: Vec<&BigInt> = roots[1..].iter().rev().collect();
        let mut inv_fft = self.compute_fft(&values, &roots);
        inv_fft.iter_mut().for_each(|inv| {
            *inv = &*inv % &self.p;
        });
        return inv_fft;
    }
}

#[cfg(test)]
mod test {
    use num_bigint::BigInt;

    use crate::fft::FFT;

    #[test]
    fn test_root_of_unity() {
        let root_of_unity = BigInt::from(9);
        let p = BigInt::from(17);
        let fft = FFT { p };
        let result = fft.build_root_of_unity(&root_of_unity);

        println!("{:?}", result)
    }
}
