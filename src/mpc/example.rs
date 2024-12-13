extern crate nalgebra as na;

use crate::mpc::linear_discrete_model::LinearDiscreteModel;

#[derive(Clone, Copy)]
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

pub fn generate_pulse_trajectory(time_steps: usize) -> na::DMatrix::<f64>
{
    let mut traj = na::DMatrix::<f64>::zeros(time_steps, 1);
    traj
        .view_range_mut(0..time_steps / 3, 0..1)
        .copy_from(&na::DMatrix::from_element(time_steps / 3, 1, 1.0));
    traj
        .view_range_mut(2 * time_steps / 3..time_steps, 0..1)
        .copy_from(&na::DMatrix::from_element(time_steps / 3, 1, 1.0));

    traj
}
