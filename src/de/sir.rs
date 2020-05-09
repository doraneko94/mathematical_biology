/*
use gnuplot::AxesCommon;
use gnuplot::*;

pub struct SirModel {
    pub beta: f64,
    pub b: f64,
    pub gamma: f64,
    pub save_t: Vec<usize>,
    pub save_s: Vec<f64>,
    pub save_i: Vec<f64>,
    pub save_r: Vec<f64>,
}

impl SirModel {
    pub fn new(beta: f64, b: f64, gamma: f64) -> Self {
        let save_t = Vec::new();
        let save_s = Vec::new();
        let save_i = Vec::new();
        let save_r = Vec::new();
        SirModel { beta, b, gamma, save_t, save_s, save_i, save_r, }
    }

    pub fn run(&mut self, time: usize, s0: f64, i0: f64) {
        let n = s0 + i0;
        let mut s = s0;
        let mut i = i0;
        self.save_t = Vec::with_capacity(time+1);
        self.save_s = Vec::with_capacity(time+1);
        self.save_i = Vec::with_capacity(time+1);
        self.save_r = Vec::with_capacity(time+1);
        self.save_t.push(0);
        self.save_s.push(s);
        self.save_i.push(i);
        self.save_r.push(0.0);

        for t in 1..time+1 {
            let (s1, i1) = self.f(s, i, n);
            s = s1;
            i = i1;
            self.save_t.push(t);
            self.save_s.push(s);
            self.save_i.push(i);
            self.save_r.push(n-s-i);
        }
    }

    fn f(&self, s0: f64, i0: f64, n: f64) -> (f64, f64) {
        let s1 = s0 - self.beta * i0 * s0 / n + self.b * (n - s0);
        let i1 = i0 * (1.0 - self.gamma - self.b) + self.beta * i0 * s0 / n;
        (s1, i1)
    }

    pub fn plot(&self, savename: &str, pos_legend: (f64, f64)) {
        let time = self.save_t.len() - 1;
        let n = self.save_s[0] + self.save_i[0];

        let mut fg = gnuplot::Figure::new();
        fg.axes2d()
            .points(self.save_t.iter(), self.save_s.iter(), &[gnuplot::PointSymbol('o'), gnuplot::Color("blue"), Caption("S (Susceptible)")])
            .points(self.save_t.iter(), self.save_i.iter(), &[gnuplot::PointSymbol('o'), gnuplot::Color("red"), Caption("I (Infected)")])
            .points(self.save_t.iter(), self.save_r.iter(), &[gnuplot::PointSymbol('o'), gnuplot::Color("green"), Caption("R (Removed)")])
            .set_x_label("Year", &[])
            .set_y_label("S, I, R", &[])
            .set_title("SIR model", &[])
            .set_legend(Coordinate::Graph(pos_legend.0), Coordinate::Graph(pos_legend.1), &[], &[])
            .set_x_range(Fix(0.0), Fix(time as f64))
            .set_y_range(Fix(0.0), Fix(n));
        fg.echo_to_file(savename);
    }

    pub fn r0(&self) {
        let r0 = self.beta / (self.b + self.gamma);
        if r0 < 1.0 {
            println!("R0 = {:.03} < 1", r0);
        } else if r0 == 1.0 {
            println!("R0 = {:.03} = 1", r0);
        } else {
            println!("R0 = {:.03} > 1", r0);
        }

        println!("");
        println!("=======================================================================");
        println!("infection     | place (period)                            | R0");
        println!("-----------------------------------------------------------------------");
        println!("smallpox      | developing countries (before eradication) |  3.0 -  5.0");
        println!("measles       | England & Wales (1956-1968)               |        13.0");
        println!("              | United States (1910-1930)                 | 12.0 - 13.0");
        println!("pertussis     | England & Wales (1942-1950)               |        17.0");
        println!("              | Maryland state (1908-1917)                |        13.0");
        println!("rubella       | England & Wales (1979)                    |         6.0");
        println!("varicella     | United States (1913-1921, 1943)           |  9.0 - 10.0");
        println!("diphteria     | United States (1910-1947)                 |  4.0 -  6.0");
        println!("scarlatina    | United States (1910-1920)                 |  5.0 -  7.0");
        println!("mumps         | United States (1912-1916, 1943)           |  4.0 -  7.0");
        println!("poliomyelitis | Netherlands(1960) & United States (1955)  |         6.0");
        println!("=======================================================================");
        println!("                                                             May (1983)")
    }
}
*/
use eom::traits::*;
use ndarray::*;

#[derive(Clone, Copy, Debug)]
pub struct Sir {
    pub beta: f64,
    pub b: f64,
    pub gamma: f64,
}

impl Default for Sir {
    fn default() -> Self {
        Sir {
            beta: 0.3,
            b: 0.2,
            gamma: 0.2,
        }
    }
}

impl ModelSpec for Sir {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        3
    }
}

impl Sir {
    pub fn new(beta: f64, b: f64, gamma: f64) -> Self {
        Sir { beta, b, gamma, }
    }
}

impl Explicit for Sir {
    fn rhs<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = f64>,
    {
        let s = v[0];
        let i = v[1];
        let r = v[2];
        let n = s + i + r;
        v[0] = s - self.beta * i * s / n + self.b * (i + r);
        v[1] = i * (1.0 - self.gamma - self.b) + self.beta * i * s / n;
        v[2] = n - v[0] - v[1];
        v
    }
}