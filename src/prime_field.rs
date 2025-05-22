use std::str::FromStr;

use num_bigint::{BigInt, Sign};

/// Represent a Prime field and implementation of operations that can be done in the prime field
pub struct PrimeField {
    p: BigInt,
}

impl PrimeField {
    pub fn add(&self, a: &BigInt, b: &BigInt) -> BigInt {
        (a + b) % &self.p
    }

    pub fn sub(&self, a: &BigInt, b: &BigInt) -> BigInt {
        (a - b) % &self.p
    }

    pub fn mul(&self, a: &BigInt, b: &BigInt) -> BigInt {
        (a * b) % &self.p
    }

    pub fn div(&self, a: &BigInt, b: &BigInt) -> BigInt {
        self.mul(a, &self.mod_inverse(&b))
    }

    /// Modular inverse using euclidean extended algorithm
    /// Modular inverse for prime always exist for a > 0
    pub fn mod_inverse(&self, a: &BigInt) -> BigInt {
        a.modinv(&self.p).unwrap()

        // if a == BigInt::ZERO {
        //     return BigInt::ZERO;
        // }

        // let mut t = BigInt::ZERO;
        // let mut new_t = BigInt::from(1);

        // let mut r = self.p;
        // let (_, mut new_r) = a.div_mod(self.p);

        // while new_r > U256::one() {
        //     let (quotient, _) = r.div_mod(new_r);
        //     t = new_t;
        //     println!("Quot {}", quotient);
        //     new_t = t - quotient * new_t;

        //     r = new_r;
        //     new_r = r - quotient * new_r;
        // }

        // let (_, t) = t.div_mod(self.p);
        // t
    }

    /// Montgomery batch inversion
    pub fn batch_inverse(&self, values: &[BigInt]) -> Vec<BigInt> {
        let n = values.len();
        let mut partials = Vec::with_capacity(n + 1);
        partials.push(BigInt::from(1));

        for value in values {
            let last = partials.last().unwrap();
            partials.push(self.mul(last, value));
        }

        let mut inv = self.mod_inverse(&partials[n]);
        let mut outputs = vec![BigInt::from(0); n];

        for i in (0..n).rev() {
            if values[i] != BigInt::from(0) {
                outputs[i] = self.mul(&partials[i], &inv);
                inv = self.mul(&inv, &values[i]);
            }
        }

        outputs
    }

    /// Function to evaluate the given polynomial at x in the prime field
    /// Polynomial is represented in array where each array is the coefficient of degree in ascending order
    pub fn evaluate_polynomial(&self, polynomial: &[BigInt], x: BigInt) -> BigInt {
        let mut x_evaluated = BigInt::from(1);
        let mut y = BigInt::ZERO;
        for p in polynomial {
            y += x_evaluated.clone() * p;

            // Evaluate X for the next iteration
            x_evaluated = (x_evaluated * x.clone()) % self.p.clone();
        }

        y % self.p.clone()
    }

    pub fn add_polynomial(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        let max_len = std::cmp::max(a.len(), b.len());
        let mut result = Vec::with_capacity(max_len);
        let zero = BigInt::ZERO;
        for i in 0..max_len {
            let ai = a.get(i).unwrap_or(&zero);
            let bi = b.get(i).unwrap_or(&zero);
            result.push((ai + bi) % self.p.clone());
        }

        result
    }

    pub fn sub_polynomial(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        let max_len = std::cmp::max(a.len(), b.len());
        let mut result = Vec::with_capacity(max_len);
        let zero = BigInt::ZERO;
        for i in 0..max_len {
            let ai = a.get(i).unwrap_or(&zero);
            let bi = b.get(i).unwrap_or(&zero);
            result.push((ai - bi) % self.p.clone());
        }

        result
    }

    pub fn mul_polynomial(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        let mut out = vec![BigInt::ZERO; a.len() + b.len() - 1];
        for i in 0..a.len() {
            for j in 0..b.len() {
                out[i + j] += (a[i].clone() + b[j].clone()) % self.p.clone();
            }
        }

        out
    }

    pub fn div_polynomial(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        let mut a = a.to_vec();
        assert!(a.len() >= b.len() && b.len() != 0);
        let mut result = Vec::new();
        let mut pos_a = a.len() - 1;
        let pos_b = b.len() - 1;
        let mut diff = (pos_a - pos_b) as isize;
        while diff >= 0 {
            //Divide the first term of dividend by highest term of divisor
            let quotient = self.div(&a[pos_a], &b[pos_b]);

            // Multiply the quotient and the divisor and eliminiate the fist element of the dividend (Because we don't that anymore)
            for i in (0..pos_b).rev() {
                a[diff as usize + i] -= &quotient * &b[i];
            }
            pos_a -= 1;
            diff -= 1;
            result.push(quotient);
        }
        result.reverse();
        result.iter().map(|r| r % &self.p).collect()
    }

    pub fn mod_polynomial(&self, a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
        self.sub_polynomial(a, &self.mul_polynomial(b, &self.div_polynomial(a, b)))
            .into_iter()
            .take(b.len() - 1)
            .collect()
    }

    /// Returns a polynomial which returns 0 at all x (vanishing polynomial)
    ///
    /// Example
    /// (x-2).(x-3).(x-4) for x in (2,3,4)
    pub fn zeropoly(&self, xs: &[BigInt]) -> Vec<BigInt> {
        // 1 here is c in the polynomial ax^2 + bx + c
        // We insert next x degree to the left side of the array
        let mut root = vec![BigInt::from(1)];
        let len = root.len();
        for x in xs {
            root.insert(0, BigInt::ZERO);
            for i in 0..len - 1 {
                let med = root[i + 1].clone() * x;
                root[i] -= med;
            }
        }
        root.iter().map(|r| r % &self.p).collect()
    }

    /// Interpolate polynomial using lagrange
    pub fn lagrange_interpolation(&self, xs: &[BigInt], ys: &[BigInt]) -> Vec<BigInt> {
        let root = self.zeropoly(xs);

        let mut numerators = Vec::with_capacity(xs.len());

        for x in xs {
            numerators.push(self.div_polynomial(&root, &[-x, BigInt::from(1)]));
        }

        //Generate denominators by evaluating numerator at each x
        let mut denominators = Vec::with_capacity(xs.len());
        for (i, x) in xs.iter().enumerate() {
            denominators.push(self.evaluate_polynomial(&numerators[i], x.clone()));
        }

        let denoms_inverse = self.batch_inverse(&denominators);

        // Generate output polynomial, which is the sum of the per-value numerator
        // polynomials rescaled to have the right y values
        let mut b = vec![BigInt::ZERO; ys.len()];

        for i in 0..xs.len() {
            let yslice = self.mul(&ys[i], &denoms_inverse[i]);
            for j in 0..ys.len() {
                if numerators[i][j] > BigInt::ZERO && ys[i] > BigInt::ZERO {
                    b[j] += &numerators[i][j] * &yslice;
                }
            }
        }

        b.iter().map(|b| b % &self.p).collect()
    }
}

#[test]
fn test_prime() {
    let big = BigInt::from_str(
        "115792089237316195423570985008687907852837564279074904382605163141518161494017",
    )
    .unwrap();

    let v = BigInt::from(31);
    let prime = PrimeField { p: v };

    let poly = [BigInt::from(4), BigInt::from(5), BigInt::from(6)];
    let m = prime.evaluate_polynomial(&poly, BigInt::from(2));

    println!("{:?}", m)
}

