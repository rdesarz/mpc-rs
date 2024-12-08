use crate::mpc::simulator::system_simulate;
use ndarray::{array, Array, Array2};
use ndarray_linalg::{Eig, Inverse};

pub fn compute_system_response() -> Vec<f64> {
    // Define parameters
    let m1 = 2.0;
    let m2 = 2.0;
    let k1 = 100.0;
    let k2 = 200.0;
    let d1 = 1.0;
    let d2 = 5.0;

    // Define the continuous-time system matrices
    let mat_ac = array![
        [0.0, 1.0, 0.0, 0.0],
        [-(k1 + k2) / m1, -(d1 + d2) / m1, k2 / m1, d2 / m1],
        [0.0, 0.0, 0.0, 1.0],
        [k2 / m2, d2 / m2, -k2 / m2, -d2 / m2]
    ];
    let mat_bc = array![[0.0], [0.0], [0.0], [1.0 / m2]];
    let mat_cc = array![[1.0, 0.0, 0.0, 0.0]];

    // Discretization constant
    let sampling = 0.05;

    // Model discretization
    let mat_i: Array2<f64> = Array::eye(mat_ac.nrows());
    let mut mat_a: Array2<f64> = mat_i - sampling * mat_ac.clone();
    mat_a = mat_a.inv().unwrap();
    let mat_b = mat_a.dot(&(sampling * mat_bc));
    let mat_c = mat_cc;

    // check the eigenvalues
    let _eigen_a = mat_ac.eig().unwrap();
    let _eigen_aid = mat_a.eig().unwrap();

    let time_sample_test = 200;

    // Compute the system's step response
    let input_test = 10.0 * Array2::ones((1, time_sample_test));
    let x0_test = Array::zeros(4);

    // // # simulate the discrete-time system
    let (y_test, _x_test) = system_simulate(&mat_a, &mat_b, &mat_c, &input_test, &x0_test);

    let system_response : Vec<f64> = y_test.row(0).iter().map(|&val| val as f64).collect();

    system_response
}