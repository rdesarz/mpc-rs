
use ndarray::{s, Array1, Array2};

pub fn system_simulate(
    mat_a: &Array2<f64>,
    mat_b: &Array2<f64>,
    mat_c: &Array2<f64>,
    mat_u: &Array2<f64>,
    x0: &Array1<f64>,
) -> (Array2<f64>, Array2<f64>) {
    let sim_time = mat_u.shape()[1];
    let n = mat_a.shape()[0];
    let r = mat_c.shape()[0];
    let mut mat_x = Array2::zeros((n, sim_time + 1));
    let mut mat_y = Array2::zeros((r, sim_time));
    for i in 0..sim_time {
        if i == 0 {
            mat_x.slice_mut(s![.., i]).assign(&x0);
            mat_y.slice_mut(s![.., i]).assign(&(mat_c.dot(x0)));
            mat_x
                .slice_mut(s![.., i + 1])
                .assign(&(mat_a.dot(x0) + mat_b.dot(&mat_u.slice(s![.., i]))));
        } else {
            mat_y
                .slice_mut(s![.., i])
                .assign(&(mat_c.dot(&mat_x.slice(s![.., i]))));

            let mat_x_slice = mat_x.slice(s![.., i]).to_owned();
            mat_x
                .slice_mut(s![.., i + 1])
                .assign(&(mat_a.dot(&mat_x_slice) + mat_b.dot(&mat_u.slice(s![.., i]))));
        }
    }

    (mat_y, mat_x)
}
