extern crate nalgebra as na;

use crate::mpc::linear_discrete_model::LinearDiscreteModel;
use crate::mpc::simulator::system_simulate;

pub struct Model {
    mat_a: na::DMatrix<f64>,
    mat_b: na::DMatrix<f64>,
    mat_c: na::DMatrix<f64>,
    sampling_dt: f64,
}

impl Model {
    pub fn new(m1: f64, m2: f64, k1: f64, k2: f64, d1: f64, d2: f64, sampling_dt: f64) -> Model {
        // Define the continuous-time system matrices
        let mat_ac = na::dmatrix![
            0.0, 1.0, 0.0, 0.0;
            -(k1 + k2) / m1, -(d1 + d2) / m1, k2 / m1, d2 / m1;
            0.0, 0.0, 0.0, 1.0;
            k2 / m2, d2 / m2, -k2 / m2, -d2 / m2
        ];
        let mat_bc = na::dmatrix![0.0; 0.0; 0.0; 1.0 / m2];
        let mat_cc = na::dmatrix![1.0, 0.0, 0.0, 0.0];

        // Model discretization
        let mat_i = na::DMatrix::<f64>::identity(mat_ac.nrows(), mat_ac.nrows());
        let mat_a = (mat_i - mat_ac.scale(sampling_dt)).try_inverse().unwrap();
        let mat_b = &mat_a * mat_bc.scale(sampling_dt);
        let mat_c = mat_cc;

        Model {
            mat_a: mat_a,
            mat_b: mat_b,
            mat_c: mat_c,
            sampling_dt: sampling_dt,
        }
    }
}

impl LinearDiscreteModel for Model {
    fn get_mat_a(&self) -> &na::DMatrix<f64> {
        &self.mat_a
    }

    fn get_mat_b(&self) -> &na::DMatrix<f64> {
        &self.mat_b
    }

    fn get_mat_c(&self) -> &na::DMatrix<f64> {
        &self.mat_c
    }

    fn get_sampling_dt(&self) -> f64 {
        self.sampling_dt
    }
}

pub fn compute_system_response(model: &impl LinearDiscreteModel, sampling_time: f64) -> Vec<f64> {
    // Compute the system's step response
    let n_samples = (sampling_time / model.get_sampling_dt()).floor() as usize;

    let input_test = na::DMatrix::from_element(1, n_samples, 10.0f64);
    let x0_test = na::DVector::<f64>::zeros(4);

    // // # simulate the discrete-time system
    let (y_test, _x_test) = system_simulate(
        &model.get_mat_a(),
        &model.get_mat_b(),
        &model.get_mat_c(),
        &input_test,
        &x0_test,
    );

    let system_response: Vec<f64> = y_test.row(0).iter().map(|&val| val as f64).collect();

    system_response
}
