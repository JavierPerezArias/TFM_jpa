mod shor;

fn miller_test(n: u32) -> bool {
    let bases = vec![2, 3, 5, 7, 11, 13, 17];
    let n1: u32 = n - 1;
    let s: u32 = n1.trailing_zeros();
    let d: u32 = n1 >> s;

    if bases.contains(&n) {
        return true;
    }
    for a in bases {
        let mut x = shor::sim::utils::mod_pow(a, d, n);
        let mut y = 1;
        for _ in 0..s {
            y = shor::sim::utils::mod_pow(x, 2, n);
            if (y == 1) && (x != 1) && (x != n1) {
                return false;
            }
            x = y;
        }
        if y != 1 {
            return false;
        }
    }

    true
}

fn is_prime(num: u32) -> bool {
    match num {
        n if n < 2 => false,
        n if n == 2 => true,
        n if (n & 1) == 0 => false,
        n => miller_test(n),
    }
}

fn base_exp(num: u32) -> Option<(u32, u32)> {
    for base in 2..format!("{:b}", num).chars().count() as u32 {
        let exp = (num as f64).powf(1.0 / (base as f64));
        if (exp - exp.round()).abs() < 1e-10 {
            return Some((exp as u32, base));
        }
    }

    None
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if !(2..=3).contains(&args.len()) {
        return;
    }
    let num: u32 = args[1].parse().unwrap();
    let t: u32 = if args.len() < 3 {
        3
    } else {
        args[2].parse().unwrap()
    };

    match num {
        n if n < 2 => println!("N={}. N must be greater than 1.", n),
        n if is_prime(n) => println!("N={} is prime.", n),
        n if (n & 1) == 0 => println!("N={} is even. p={}, q={}.", n, 2, n / 2),
        n => {
            if let Some((base, exp)) = base_exp(n) {
                println!("N={}. N is exponentiation. {}^{}.", n, base, exp)
            } else {
                shor::shor(num, t)
            }
        }
    }
}
