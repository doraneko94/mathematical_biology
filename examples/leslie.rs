use mathematical_biology::leslie::*;
use numerical::difference::*;
use ndarray::*;
use eom::*;

fn main() {
    let eom = Leslie::default();
    let mut teo = Diff::from(eom.clone());
    let ts = adaptor::time_series(arr1(&[24.0, 4.0, 1.0]), &mut teo);
    let end_time = 10;
    
    println!("*** Default ***");
    println!("L = [[   0,   9,  12],");
    println!("     [ 1/3,   0,   0],");
    println!("     [   0, 1/2,   0]]");
    println!("t =  0: x0 =    24, x1 =     4, x2 =     1");
    for (i, v) in ts.take(end_time).enumerate() {
        println!("t = {:>2}: x0 = {:>5.0}, x1 = {:>5.0}, x2 = {:>5.0}",
        i+1, v[0], v[1], v[2]);
    }
    println!("** stats **");
    println!("lambda1 = {:>7.5}", eom.lambda1);
    println!("V1      = {:>7.5}", eom.v1);
    println!("R0      = {:>7.5}", eom.r0);
    println!("");

    let eom = Leslie::new(&[[0.0000, 0.0000, 0.0000, 0.0000, 127.00, 4.0000, 80.000],
                            [0.6747, 0.7370, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
                            [0.0000, 0.0486, 0.6610, 0.0000, 0.0000, 0.0000, 0.0000],
                            [0.0000, 0.0000, 0.0147, 0.6907, 0.0000, 0.0000, 0.0000],
                            [0.0000, 0.0000, 0.0000, 0.0518, 0.0000, 0.0000, 0.0000],
                            [0.0000, 0.0000, 0.0000, 0.0000, 0.8091, 0.0000, 0.0000],
                            [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.8091, 0.8089]]);
    let mut teo = Diff::from(eom.clone());
    let s = [20651.0, 66975.0, 11460.0, 662.0, 36.0, 31.0, 185.0];
    let total0: f64 = s.iter().sum();
    let ts = adaptor::time_series(arr1(&s), &mut teo);
    let end_time = 25;

    println!("*** Caretta caretta [Crouse et al. 1987] ***");
    println!("L = [[0.0000, 0.0000, 0.0000, 0.0000, 127.00, 4.0000, 80.000],");
    println!("     [0.6747, 0.7370, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],");
    println!("     [0.0000, 0.0486, 0.6610, 0.0000, 0.0000, 0.0000, 0.0000],");
    println!("     [0.0000, 0.0000, 0.0147, 0.6907, 0.0000, 0.0000, 0.0000],");
    println!("     [0.0000, 0.0000, 0.0000, 0.0518, 0.0000, 0.0000, 0.0000],");
    println!("     [0.0000, 0.0000, 0.0000, 0.0000, 0.8091, 0.0000, 0.0000],");
    println!("     [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.8091, 0.8089]]");
    println!("t =  0: x0 = 20651, x1 = 66975, x2 = 11460, x3 =   662, x4 =    36, x5 =    31, x6 =   185, total = {:>10.4}",
             total0);
    let mut total25 = 0.0;
    for (i, v) in ts.take(end_time).enumerate() {
        println!("t = {:>2}: x0 = {:>5.0}, x1 = {:>5.0}, x2 = {:>5.0}, x3 = {:>5.0}, x4 = {:>5.0}, x5 = {:>5.0}, x6 = {:>5.0}, total = {:>10.4}",
        i+1, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v.sum());
        if i == 24 { total25 = v.sum() }
    }
    println!("** stats **");
    println!("lambda1 = {:>7.5}", eom.lambda1);
    println!("V1      = {:>7.5}", eom.v1);
    println!("R0      = {:>7.5}", eom.r0);
    println!("** result **");
    let lambda1_25 = eom.lambda1.powi(25);
    println!("(lambda1)^25 = {:>7.5}", lambda1_25);
    println!("total0 * (lambda1)^25 = {:>7.5} ~ {:>7.5} = total25", lambda1_25 * total0, total25);
}