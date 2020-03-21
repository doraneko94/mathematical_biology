use gnuplot::AxesCommon;
use gnuplot::*;

fn f(x: f64, r: f64) -> f64 {
    r * x * ( 1.0 - x )
}

fn main() {
    let mut x = Vec::new();
    let mut y = Vec::new();

    let mut r = 3.0;
    while r < 4.0 {
        let mut xt = 0.5;
        for _ in 0..400 {
            xt = f(xt, r);
        }
        let mut s = 0.0;
        for _ in 0..500 {
            xt = f(xt, r);
            s += ( r - 2.0 * r * xt ).abs().log(std::f64::consts::E);
        }
        x.push(r);
        y.push(s / 500.0);
        r += 0.001;
    }
    
    let mut fg = gnuplot::Figure::new();
    fg.axes2d()
        .lines(x.iter(), y.iter(), &[gnuplot::Color("blue")])
        .set_x_label("r", &[])
        .set_y_label("Liapunov Index", &[])
        .set_title("Liapunov", &[])
        .set_x_range(Fix(3.0), Fix(4.0))
        .set_y_range(Fix(-4.0), Fix(1.0));
    fg.echo_to_file("liapunov.plt");
}