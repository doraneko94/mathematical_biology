use ndarray::*;
use ndarray_linalg::*;
use numerical::difference::Difference;
use eom::traits::*;

#[derive(Clone, Debug)]
pub struct Leslie {
    pub l: Array2<f64>,
    pub lambda1: f64,
    pub v1: Array1<f64>,
    pub r0: f64,
}

impl Default for Leslie {
    fn default() -> Self {
        Self::new(&[[0.0, 9.0, 12.0],
                    [1.0/3.0, 0.0, 0.0],
                    [0.0, 1.0/2.0, 0.0]])
    }
}

impl Leslie {
    pub fn new<V: FixedInitializer<Elem = f64>>(l: &[V]) -> Self
    where
        V: Clone
    {
        let l = arr2(l);
        let d = l.shape();
        if d[0] != d[1] {
            panic!("Row and Column should have equal sizes!");
        }
        let (lambda1, v1, r0) = calc_eig(&l);
        Self { l, lambda1, v1, r0, }
    }
}

impl ModelSpec for Leslie {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        self.v1.shape()[0]
    }
}

impl Difference for Leslie {
    fn recr<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = Self::Scalar>
    {
        let new_v = self.l.dot(v);

        for i in 0..new_v.shape()[0] {
            v[i] = new_v[i];
        }
        v
    }
}

fn calc_eig(l: &Array2<f64>) -> (f64, Array1<f64>, f64) {
    let (e, vecs) = l.clone().eig().unwrap();
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
    let v_scale: Vec<f64>  = v.iter().map(|e| e / v_sum).collect();

    let mut r0 = l[[0, 0]];
    let mut s = 1.0;
    for i in 1..n {
        s *= l[[i, i-1]];
        r0 += l[[0, i]] * s;
    }
    (e[index].re, arr1(&v_scale), r0)
}