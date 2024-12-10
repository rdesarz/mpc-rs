extern crate nalgebra as na;

use crate::mpc::simulator::system_simulate;

pub fn compute_system_response() -> Vec<f64> {
    // Define parameters
    let m1 = 2.0;
    let m2 = 2.0;
    let k1 = 100.0;
    let k2 = 200.0;
    let d1 = 1.0;
    let d2 = 5.0;
    // Define the continuous-time system matrices
    let mat_ac = na::dmatrix![
        0.0, 1.0, 0.0, 0.0;
        -(k1 + k2) / m1, -(d1 + d2) / m1, k2 / m1, d2 / m1;
        0.0, 0.0, 0.0, 1.0;
        k2 / m2, d2 / m2, -k2 / m2, -d2 / m2
    ];
    let mat_bc = na::dmatrix![0.0; 0.0; 0.0; 1.0 / m2];
    let mat_cc = na::dmatrix![1.0, 0.0, 0.0, 0.0];

    // Discretization constant
    let sampling = 0.05f64;

    // Model discretization
    let mat_i = na::DMatrix::<f64>::identity(mat_ac.nrows(), mat_ac.nrows());
    let mat_a = (mat_i - mat_ac.scale(sampling)).try_inverse().unwrap();
    let mat_b = &mat_a * mat_bc.scale(sampling);
    let mat_c = mat_cc;

    let time_sample_test = 200;

    // Compute the system's step response
    let input_test = na::DMatrix::from_element(1, time_sample_test, 10.0f64);
    let x0_test = na::DVector::<f64>::zeros(4);

    // // # simulate the discrete-time system
    let (y_test, _x_test) = system_simulate(&mat_a, &mat_b, &mat_c, &input_test, &x0_test);

    let system_response : Vec<f64> = y_test.row(0).iter().map(|&val| val as f64).collect();

    system_response
}