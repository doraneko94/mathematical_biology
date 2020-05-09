use eom::traits::*;
use ndarray::*;

#[derive(Clone, Copy, Debug)]
pub struct KermackMcKendrick {
    pub beta: f64,
    pub b: f64,
    pub gamma: f64,
}

impl Default for KermackMcKendrick {
    fn default() -> Self {
        Self {
            beta: 0.3,
            b: 0.2,
            gamma: 0.2,
        }
    }
}

impl ModelSpec for KermackMcKendrick {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        3
    }
}

impl KermackMcKendrick {
    pub fn new(beta: f64, b: f64, gamma: f64) -> Self {
        Self { beta, b, gamma, }
    }
}

impl Explicit for KermackMcKendrick {
    fn rhs<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = f64>,
    {
        let s = v[0];
        let i = v[1];
        let ds = self.beta * s * i;
        let di = self.gamma * i;
        v[0] = - ds;
        v[1] = - di + ds;
        v[2] = di;

        v
    }
}