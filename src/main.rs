use num_bigint::BigInt;

mod fft;
mod fri;
mod prime_field;

fn main() {
    println!("Hello, world!");
}

fn find_root_of_unity(prime: BigInt) {
    let mut p = prime - BigInt::from(1);

    //Factorise p using 2
    let mut k = 0; // 2^k
    while &p % BigInt::from(2) == BigInt::from(0) {
        p = p / BigInt::from(2);

        k += 1;
    }
    
    //Find a primitive generator g
    let mut g = BigInt::from(2);
    // while g.pow(p) != BigInt::from(1) {
    //     g = g + BigInt::from(1);
    // }
}
