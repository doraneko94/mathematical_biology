use ndarray::*;
use ndarray_linalg::*;
use num_complex::Complex64 as c64;

fn main() {
    let l = arr2(&[[0.0000, 0.0000, 0.0000, 0.0000, 127.00, 4.0000, 80.000],
                   [0.6747, 0.7370, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
                   [0.0000, 0.0486, 0.6610, 0.0000, 0.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, 0.0147, 0.6907, 0.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, 0.0000, 0.0518, 0.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, 0.0000, 0.0000, 0.8091, 0.0000, 0.0000],
                   [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.8091, 0.8089]]);
    let (e, vecs) = l.clone().eig().unwrap();

    let (_l1, v1) = l1v1(&e, &vecs);
    println!("{:?}", v1);
}

fn l1v1(e: &Array1<c64>, vecs: &Array2<c64>) -> (f64, Vec<f64>) {
    let n = e.raw_dim()[0];
    let mut index = 0;
    let mut m = 0.0;
    for (i, l) in e.iter().enumerate() {
        if m < l.abs() {
            m = l.abs();
            index = i;
        }
    }
    if e[index].im != 0.0 {
        panic!("Dominant eigenvalue shouldn't be a complex number!");
    }
    let v: Vec<f64> = vecs.slice(s![0..n, index]).iter().map(|c| c.re).collect();
    let v_sum: f64 = v.iter().sum();
    let v = v.iter().map(|e| e / v_sum).collect();
    (e[index].re, v)
}