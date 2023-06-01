use num::{Complex, Zero};
use std::rc::Rc;

pub mod utils;

type BinaryChars = Vec<char>;
type Qubit = Vec<Complex<f64>>;
type Gate = Vec<Qubit>;

pub struct State {
    number_of_bit: u32,
    index: usize,
    amp: Complex<f64>,
    prob: f64,
}
impl State {
    pub fn to_binary_chars(&self, qb: &[u32]) -> BinaryChars {
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

fn tensor_with(list: &[Rc<Gate>]) -> Gate {
    let mut g = list[0].to_vec();

    for i in list.iter().skip(1) {
        let mut t = vec![];
        let n = Rc::clone(i);

        for j in 0..g.len() {
            for k in 0..n.len() {
                let mut v = vec![];

                for l in 0..g[j].len() {
                    for m in 0..n[k].len() {
                        v.push(g[j][l] * n[k][m]);
                    }
                }
                t.push(v);
            }
        }
        g = t;
    }

    g
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

fn cmodexp2(nob: u32, a: u32, j: u32, n: u32, control: u32, target: &[u32]) -> Gate {
    let r1len = target.len() as u32;
    let r0len = nob - r1len;
    let a2jmodn = utils::mod_pow(a, j.pow(2), n);
    let mut index = vec![];

    for i in 0..(2_usize.pow(nob)) {
        let bits: BinaryChars = format!("{:>0n$b}", i, n = nob as usize).chars().collect();

        if bits[control as usize] == '0' {
            index.push(utils::to_decimal(&bits) as usize);
            continue;
        }
        let k = utils::to_decimal(&bits[r0len as usize..bits.len()].to_vec());
        if k > n - 1 {
            index.push(utils::to_decimal(&bits) as usize);
            continue;
        }
        let mut a2jkmodns: BinaryChars =
            format!("{:>0n$b}", ((a2jmodn * k) % n) as usize, n = r1len as usize)
                .chars()
                .collect();
        let mut r0bits = bits[0..r0len as usize].to_vec();
        r0bits.append(&mut a2jkmodns);
        index.push(utils::to_decimal(&r0bits) as usize);
    }
    let id = id_with(nob);
    let mut g = vec![vec![]; id.len()];
    for (i, ii) in index.iter().enumerate() {
        g[*ii] = id[i].to_vec();
    }

    g
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

pub struct Q {
    qb: Qubit,
}
impl Q {
    pub fn new() -> Q {
        Q { qb: vec![] }
    }

    pub fn zero_with(&mut self, n: u32) -> Vec<u32> {
        let mut list = vec![];

        for _ in 0..n {
            let qb = vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
            if self.qb.is_empty() {
                self.qb = qb;
                list.push(0);
                continue;
            }
            let mut v = vec![];
            for w in &self.qb {
                for j in &qb {
                    v.push(w * j);
                }
            }
            self.qb = v;
            list.push((self.qb.len() as f64).log2() as u32 - 1);
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

    pub fn cmodexp2(&mut self, a: u32, n: u32, r0: &[u32], r1: &[u32]) {
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

    pub fn x(&mut self, qb: &[u32]) {
        self.apply_with(
            vec![
                vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            ],
            qb,
        )
    }

    pub fn cr(&mut self, theta: f64, control: u32, target: u32) {
        self.apply(cr(
            theta,
            (self.qb.len() as f64).log2() as u32,
            control,
            target,
        ))
    }

    pub fn h(&mut self, qb: &[u32]) {
        let hdm = Complex::new(1.0 / std::f64::consts::SQRT_2, 0.0);

        self.apply_with(vec![vec![hdm, hdm], vec![hdm, -1.0 * hdm]], qb)
    }

    pub fn iqft(&mut self, qb: &[u32]) {
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

    pub fn state(&self) -> Vec<State> {
        let mut list = vec![];

        for i in 0..self.qb.len() {
            let rqb = utils::round(self.qb[i]);

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
