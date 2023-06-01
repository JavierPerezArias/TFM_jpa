use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand_pcg::Pcg32;

pub mod sim;

const SEED: u64 = 1;

fn parse_float(bin: &[char]) -> f64 {
    let mut f = 0.0;

    for (i, b) in bin.iter().enumerate() {
        if !(*b == '0') {
            f += 0.5_f64.powf((i + 1) as f64);
        }
    }

    f
}

fn continued_fraction(f: f64) -> Vec<u32> {
    let mut list = vec![];
    let mut r = f;

    loop {
        let t = r.trunc();

        list.push(t as u32);
        let diff = r - t;
        if diff < 1e-3 {
            break;
        }
        r = 1.0 / diff;
    }

    list
}

fn convergent(cf: &[u32]) -> (u32, u32) {
    let len = cf.len();

    if len == 1 {
        return (cf[0], 1);
    }
    let mut s = 1;
    let mut r = cf[len - 1];
    for i in 2..len {
        let tmp = s;
        s = r;
        r = cf[len - i] * r + tmp;
    }
    s += cf[0] * r;

    (s, r)
}

fn find_order(a: u32, n: u32, bin: &[char]) -> (u32, u32, bool) {
    if bin.is_empty() {
        return (0, 1, false);
    }
    let fv = parse_float(bin);
    let cf = continued_fraction(fv);

    for i in 0..cf.len() {
        let (s, r) = convergent(&cf[0..(i + 1)]);

        if sim::utils::mod_pow(a, r, n) == 1 {
            return (s, r, true);
        }
    }
    let (s, r) = convergent(&cf);

    (s, r, false)
}

fn is_trivial(n: u32, factor: &[u32]) -> bool {
    for (_, p) in factor.iter().enumerate() {
        if 1 < *p && *p < n && n % p == 0 {
            return false;
        }
    }

    true
}

pub fn shor(num: u32, t: u32) {
    let mut rng = Pcg32::seed_from_u64(SEED);
    let dist = Uniform::from(2..num);
    let mut used = vec![];

    'LOOP: loop {
        let mut a = dist.sample(&mut rng);

        while used.contains(&a) {
            let r = sim::utils::gcd(a, num);
            if r != 1 {
                println!("p={}, q={} [!]", r, a);
                return;
            }
            a = dist.sample(&mut rng);
        }
        used.push(a);
        println!("N: {} (a: {}, t: {})", num, a, t);
        let mut qsim = sim::Q::new();
        let r0 = qsim.zero_with(t);
        let r1 = qsim.zero_with(format!("{:b}", num).chars().count() as u32);
        qsim.x(&[r1[r1.len() - 1]]);
        qsim.h(&r0);
        qsim.cmodexp2(a, num, &r0, &r1);
        qsim.iqft(&r0);
        for state in qsim.state().iter() {
            let m0 = state.to_binary_chars(&r0);
            let (s, r, ok) = find_order(a, num, &m0);

            if !ok || (r & 1) == 1 {
                println!("{}; s/r={:>2}/{:>2};", state, s, r);
                continue;
            }
            let p0 = sim::utils::gcd(a.pow(r / 2) - 1, num);
            let p1 = sim::utils::gcd(a.pow(r / 2) + 1, num);
            if is_trivial(num, &[p0, p1]) {
                println!("{}; s/r={:>2}/{:>2}; [~]", state, s, r);
                continue;
            }
            println!("{}; s/r={:>2}/{:>2}; p={}, q={} [!]", state, s, r, p0, p1);
            break 'LOOP;
        }
        println!();
    }
}

