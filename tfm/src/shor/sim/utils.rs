use num::Complex;

pub fn gcd(mut a: u32, mut n: u32) -> u32 {
    while n != 0 {
        let tmp = n;
        n = a % n;
        a = tmp;
    }

    a
}

pub fn mod_pow(mut base: u32, mut exp: u32, m: u32) -> u32 {
    let mut res = 1;
    base = base % m;

    while exp > 0 {
        if (exp & 1) == 1 {
            res = (res * base) % m;
        }
        base = (base * base) % m;
        exp >>= 1;
    }

    res
}

pub fn to_decimal(v: &[char]) -> u32 {
    let s: String = v.iter().collect();

    return u32::from_str_radix(&s, 2).unwrap();
}

pub fn round(c: Complex<f64>) -> Complex<f64> {
    let mut round = c;

    if c.re.abs() < 1e-13 {
        round.re = 0.0;
    }
    if c.im.abs() < 1e-13 {
        round.im = 0.0;
    }

    round
}
