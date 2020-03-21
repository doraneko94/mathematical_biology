use mathematical_biology::sir::*;

fn main() {
    let mut model = SirModel::new(0.3, 0.2, 0.2);
    model.run(25, 70.0, 30.0);
    model.plot("sir.plt", (1.0, 0.7));
    model.r0();
}