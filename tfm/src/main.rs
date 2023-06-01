use num::Complex;
use num::Zero;
use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand_pcg::Pcg32;
use std::rc::Rc;

const SEED: u64 = 1;

type BinaryChars = Vec<char>;
type Qubit = Vec<Complex<f64>>;
type Gate = Vec<Qubit>;

struct State {
    number_of_bit: u32,
    index: usize,
    amp: Complex<f64>,
    prob: f64,
}
impl State {
    fn to_binary_chars(&self, qb: &[u32]) -> BinaryChars {
        let v: BinaryChars = format!("{:>0n$b}", self.index, n = self.number_of_bit as usize)
            .chars()
            .collect();

        let mut bin = vec![];
        for i in qb {
            bin.push(v[*i as usize]);
        }

        bin
    }
}
impl std::fmt::Display for State {
    fn fmt(&self, dest: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            dest,
            "[{:>0n$b}]({:>+.4} {:>+.4}): {:>.4}",
            self.index,
            self.amp.re,
            self.amp.im,
            self.prob,
            n = self.number_of_bit as usize,
        )
    }
}

fn tensor(m: Gate, n: Rc<Gate>) -> Gate {
    let mut g = vec![];

    for i in 0..m.len() {
        for k in 0..n.len() {
            let mut v = vec![];

            for j in 0..m[i].len() {
                for l in 0..n[k].len() {
                    v.push(m[i][j] * n[k][l]);
                }
            }
            g.push(v);
        }
    }

    g
}

fn tensor_with(list: &[Rc<Gate>]) -> Gate {
    let mut g = list[0].to_vec();

    for i in list.iter().skip(1) {
        g = tensor(g, Rc::clone(i));
    }

    g
}

fn mod_pow(mut base: u32, mut exp: u32, m: u32) -> u32 {
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

fn to_decimal(v: &[char]) -> u32 {
    let s: String = v.iter().collect();

    return u32::from_str_radix(&s, 2).unwrap();
}

fn id_with(nob: u32) -> Gate {
    let mut mat = vec![];

    for i in 0..1 << nob {
        let mut v = vec![];

        for j in 0..1 << nob {
            if i == j {
                v.push(Complex::new(1.0, 0.0));
                continue;
            }
            v.push(Complex::new(0.0, 0.0));
        }
        mat.push(v);
    }

    mat
}

fn cmodexp2(nob: u32, a: u32, j: u32, n: u32, control: u32, target: &[u32]) -> Gate {
    let r1len = target.len() as u32;
    let r0len = nob - r1len;
    let a2jmodn = mod_pow(a, j.pow(2), n);
    let mut index = vec![];

    for i in 0..(2_usize.pow(nob)) {
        let bits: BinaryChars = format!("{:>0n$b}", i, n = nob as usize).chars().collect();

        if bits[control as usize] == '0' {
            index.push(to_decimal(&bits) as usize);
            continue;
        }
        // let r1bits = take(&bits, r0len as usize, bits.len());
        let k = to_decimal(&bits[r0len as usize..bits.len()].to_vec());
        if k > n - 1 {
            index.push(to_decimal(&bits) as usize);
            continue;
        }
        let mut a2jkmodns: BinaryChars =
            format!("{:>0n$b}", ((a2jmodn * k) % n) as usize, n = r1len as usize)
                .chars()
                .collect();
        // let mut r0bits = take(&bits, 0, r0len as usize);
        let mut r0bits = bits[0..r0len as usize].to_vec();
        r0bits.append(&mut a2jkmodns);
        index.push(to_decimal(&r0bits) as usize);
    }
    let id = id_with(nob);
    let mut g = vec![vec![]; id.len()];
    for (i, ii) in index.iter().enumerate() {
        g[*ii] = id[i].to_vec();
    }

    g
}

fn gate_list(nob: u32, g: Gate, qb: &[u32]) -> Vec<Rc<Gate>> {
    let mut list = vec![];
    let id = Rc::new(id_with(1));
    let rg = Rc::new(g);

    'LOOP: for i in 0..nob {
        for j in qb {
            if i == *j {
                list.push(Rc::clone(&rg));
                continue 'LOOP;
            }
        }
        list.push(Rc::clone(&id));
    }

    list
}

fn cr(theta: f64, nob: u32, control: u32, target: u32) -> Gate {
    let mut g = id_with(nob);
    let e = Complex::new(0.0, theta).exp();

    for (i, v) in g.iter_mut().enumerate() {
        let bits: BinaryChars = format!("{:>0n$b}", i, n = nob as usize).chars().collect();

        if bits[control as usize] == '1' && bits[target as usize] == '1' {
            v[i] = e * v[i];
        }
    }

    g
}

fn round(c: Complex<f64>) -> Complex<f64> {
    let mut round = c;

    if c.re.abs() < 1e-13 {
        round.re = 0.0;
    }
    if c.im.abs() < 1e-13 {
        round.im = 0.0;
    }

    round
}

struct Q {
    qb: Qubit,
}
impl Q {
    fn new() -> Q {
        Q { qb: vec![] }
    }

    fn tensor(&mut self, qb: Qubit) {
        let mut v = vec![];

        for w in &self.qb {
            for j in &qb {
                v.push(w * j);
            }
        }
        self.qb = v
    }

    fn add(&mut self, qb: Qubit) -> u32 {
        if self.qb.is_empty() {
            self.qb = qb;
            return 0;
        }
        self.tensor(qb);

        (self.qb.len() as f64).log2() as u32 - 1
    }

    fn zero_with(&mut self, n: u32) -> Vec<u32> {
        let mut list = vec![];

        for _ in 0..n {
            list.push(self.add(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]));
        }

        list
    }

    fn apply(&mut self, g: Gate) {
        let mut v = vec![];

        for i in 0..g.len() {
            let mut e = Complex::new(0.0, 0.0);

            for j in 0..g[i].len() {
                e += g[i][j] * self.qb[j];
            }
            v.push(e);
        }
        self.qb = v
    }

    fn apply_with(&mut self, g: Gate, qb: &[u32]) {
        self.apply(tensor_with(&gate_list(
            (self.qb.len() as f64).log2() as u32,
            g,
            qb,
        )))
    }

    fn cmodexp2(&mut self, a: u32, n: u32, r0: &[u32], r1: &[u32]) {
        for (i, c) in r0.iter().enumerate() {
            self.apply(cmodexp2(
                (self.qb.len() as f64).log2() as u32,
                a,
                i as u32,
                n,
                *c,
                r1,
            ))
        }
    }

    fn x(&mut self, qb: &[u32]) {
        self.apply_with(
            vec![
                vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            ],
            qb,
        )
    }

    fn cr(&mut self, theta: f64, control: u32, target: u32) {
        self.apply(cr(
            theta,
            (self.qb.len() as f64).log2() as u32,
            control,
            target,
        ))
    }

    fn h(&mut self, qb: &[u32]) {
        let hdm = Complex::new(1.0 / std::f64::consts::SQRT_2, 0.0);

        self.apply_with(vec![vec![hdm, hdm], vec![hdm, -1.0 * hdm]], qb)
    }

    fn iqft(&mut self, qb: &[u32]) {
        let len = qb.len();

        for i in (0..len).rev() {
            let mut k = (len - i) as i32;

            for j in ((i + 1)..len).rev() {
                let theta = -2.0 * std::f64::consts::PI / (2.0_f64.powf(k as f64));

                self.cr(theta, qb[j], qb[i]);
                k -= 1;
            }
            self.h(&[qb[i]]);
        }
    }

    fn state(&self) -> Vec<State> {
        let mut list = vec![];

        for i in 0..self.qb.len() {
            let rqb = round(self.qb[i]);

            if !rqb.is_zero() {
                list.push(State {
                    number_of_bit: (self.qb.len() as f64).log2() as u32,
                    index: i,
                    amp: rqb,
                    prob: rqb.norm().powf(2.0),
                });
            }
        }

        list
    }
}

fn gcd(mut a: u32, mut n: u32) -> u32 {
    while n != 0 {
        let tmp = n;
        n = a % n;
        a = tmp;
    }

    a
}

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

        if mod_pow(a, r, n) == 1 {
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

fn shor(num: u32, t: u32) {
    let mut rng = Pcg32::seed_from_u64(SEED);
    let dist = Uniform::from(2..num);
    let mut used = vec![];

    'LOOP: loop {
        let mut a = dist.sample(&mut rng);

        while used.contains(&a) {
            let r = gcd(a, num);
            if r != 1 {
                println!("p={}, q={} [!]", r, a);
                return;
            }
            a = dist.sample(&mut rng);
        }
        used.push(a);
        println!("N: {} (a: {}, t: {})", num, a, t);
        let mut qsim = Q::new();
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
            let p0 = gcd(a.pow(r / 2) - 1, num);
            let p1 = gcd(a.pow(r / 2) + 1, num);
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

fn miller_test(n: u32) -> bool {
    let bases = vec![2, 3, 5, 7, 11, 13, 17];
    let n1: u32 = n - 1;
    let s: u32 = n1.trailing_zeros();
    let d: u32 = n1 >> s;

    if bases.contains(&n) {
        return true;
    }
    for a in bases {
        let mut x = mod_pow(a, d, n);
        let mut y = 1;
        for _ in 0..s {
            y = mod_pow(x, 2, n);
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
                shor(num, t)
            }
        }
    }
}
