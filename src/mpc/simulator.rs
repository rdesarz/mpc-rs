
use nalgebra as na;

pub fn system_simulate(
    mat_a: &na::DMatrix<f64>,
    mat_b: &na::DMatrix<f64>,
    mat_c: &na::DMatrix<f64>,
    mat_u: &na::DMatrix<f64>,
    x0: &na::DVector<f64>,
) -> (na::DMatrix<f64>, na::DMatrix<f64>) {
    let sim_time = mat_u.ncols();
    let n = mat_a.nrows();
    let r = mat_c.nrows();
    let mut mat_x = na::DMatrix::<f64>::zeros(n, sim_time + 1);
    let mut mat_y = na::DMatrix::<f64>::zeros(r, sim_time);
    for i in 0..sim_time {
        if i == 0 {
            mat_x.column_mut(i).copy_from(&x0);
            mat_y.column_mut(i).copy_from(&(mat_c * x0));
            mat_x
                .column_mut(i + 1)
                .copy_from(&(mat_a * x0 + mat_b * mat_u.column(i)));
        } else {
            mat_y
                .column_mut(i)
                .copy_from(&(mat_c * mat_x.column(i)));

            let mat_x_slice = mat_x.column(i).into_owned();
            mat_x
                .column_mut(i + 1)
                .copy_from(&(mat_a * mat_x_slice + mat_b * mat_u.column(i)));
        }
    }

    (mat_y, mat_x)
}
