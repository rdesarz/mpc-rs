use ndarray::{array, Array, Array1, Array2};
use ndarray_linalg::{Eig, Inverse};
use plotters::prelude::*;

mod simulator {
    use ndarray::{s, Array1, Array2};

    #[derive(Debug, Copy, Clone)]
    pub struct TwoWheeledSystem {
        pub x: f32,
        pub y: f32,
        pub yaw: f32,
        pub v: f32,
    }

    pub struct Command {
        pub u1: f32,
        pub u2: f32,
    }

    pub fn apply_differential_model(system: &mut TwoWheeledSystem, cmd: &Command, dt: f32) {
        let dx = system.yaw.cos() * system.v;
        let dy = system.yaw.sin() * system.v;
        let dv = cmd.u1;
        let d_yaw = system.v / 0.25 * cmd.u2.sin();

        system.x += dt * dx;
        system.y += dt * dy;
        system.yaw += dt * d_yaw;
        system.v += dt * dv;
    }

    pub fn system_simulate(
        A: &Array2<f64>,
        B: &Array2<f64>,
        C: &Array2<f64>,
        U: &Array2<f64>,
        x0: &Array1<f64>,
    ) -> (Array2<f64>, Array2<f64>) {
        let sim_time = U.shape()[1];
        let n = A.shape()[0];
        let r = C.shape()[0];
        let mut X = Array2::zeros((n, sim_time + 1));
        let mut Y = Array2::zeros((r, sim_time));
        for i in 0..sim_time {
            if i == 0 {
                X.slice_mut(s![.., i]).assign(&x0);
                Y.slice_mut(s![.., i]).assign(&(C.dot(x0)));
                X.slice_mut(s![.., i + 1])
                    .assign(&(A.dot(x0) + B.dot(&U.slice(s![.., i]))));
            } else {
                Y.slice_mut(s![.., i]).assign(&(C.dot(&X.slice(s![.., i]))));

                let X_slice = X.slice(s![.., i]).to_owned();
                X.slice_mut(s![.., i + 1])
                    .assign(&(A.dot(&X_slice) + B.dot(&U.slice(s![.., i]))));
            }
        }

        (Y, X)
    }
}

use simulator::system_simulate;
use std::env;
use debug_plotter::plot;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env::set_var("RUST_BACKTRACE", "1");

    // Define parameters
    let m1 = 2.0;
    let m2 = 2.0;
    let k1 = 100.0;
    let k2 = 200.0;
    let d1 = 1.0;
    let d2 = 5.0;
    // Define the continuous-time system matrices
    let Ac = array![
        [0.0, 1.0, 0.0, 0.0],
        [-(k1 + k2) / m1, -(d1 + d2) / m1, k2 / m1, d2 / m1],
        [0.0, 0.0, 0.0, 1.0],
        [k2 / m2, d2 / m2, -k2 / m2, -d2 / m2]
    ];
    let Bc = array![[0.0], [0.0], [0.0], [1.0 / m2]];
    let Cc = array![[1.0, 0.0, 0.0, 0.0]];

    let r = 1;
    let m = 1; // number of inputs and outputs
    let n = 4; // state dimension

    // Discretization constant
    let sampling = 0.05;

    // Model discretization
    let I: Array2<f64> = Array::eye(Ac.shape()[0]);
    let mut A: Array2<f64> = I - sampling * Ac.clone();
    A = A.inv()?;
    let B = A.dot(&(sampling * Bc));
    let C = Cc;

    // check the eigenvalues
    let eigen_A = Ac.eig()?;
    let eigen_Aid = A.eig()?;

    let time_sample_test = 200;

    // Compute the system's step response
    let input_test = 10.0 * Array2::ones((1, time_sample_test));
    let x0_test = Array1::zeros(4);

    // // # simulate the discrete-time system
    let (Y_test, X_test) = system_simulate(&A, &B, &C, &input_test, &x0_test);

    debug_plotter::plot!(Y_test where caption = "My Plot");

    Ok(())
}
