/*
use mathematical_biology::sir::*;

fn main() {
    let mut model = SirModel::new(0.3, 0.2, 0.2);
    model.run(25, 70.0, 30.0);
    model.plot("sir.plt", (1.0, 0.7));
    model.r0();
}
*/

use mathematical_biology::*;
use ndarray::arr1;
use eom::*;
use eom::traits::*;

fn main() {
    let sir = infection::sir::KermackMcKendrick::default();
    let mut teo = explicit::Euler::new(sir, 0.1);
    let ts = adaptor::time_series(arr1(&[0.7, 0.3, 0.0]), &mut teo);
    let end_time = 50;
    println!("time,s,i,r");
    for v in ts.take(end_time) {
        println!("{}, {}, {}", v[0], v[1], v[2]);
    }
}