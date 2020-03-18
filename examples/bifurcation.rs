use gnuplot::AxesCommon;
use gnuplot::*;

fn logi_eq(x: f64, r: f64) -> f64 {
    r * x * ( 1.0 - x )
}

fn run(save_name: &str, r_s: f64, r_e: f64, dr: f64, y_min: f64, y_max: f64) {
    let mut x = Vec::new();
    let mut y = Vec::new();

    let mut r = r_s;
    while r <= r_e {
        let mut xt = 0.5;
        for _ in 0..1000 {
            xt = logi_eq(xt, r);
        }
        for _ in 0..1000 {
            xt = logi_eq(xt, r);
            x.push(r);
            y.push(xt);
        }
        r += dr;
    }

    let mut fg = gnuplot::Figure::new();
    fg.axes2d()
        .points(x.iter(), y.iter(), &[gnuplot::PointSymbol('.'), gnuplot::Color("blue")])
        .set_x_label("r", &[])
        .set_y_label("Stable Cycle", &[])
        .set_title("test", &[])
        .set_x_range(Fix(r_s), Fix(r_e))
        .set_y_range(Fix(y_min), Fix(y_max));
    fg.echo_to_file(save_name);
}

fn main() {
    run("bifurcation_280-400.plt", 2.8, 4.0, 0.001, 0.0, 1.0);
    run("bifurcation_384-386.plt", 3.84, 3.86, 0.00001, 0.42, 0.58);
}