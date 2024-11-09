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

mod controller {
    use ndarray::Array2;

    pub struct Controller {
        A: Array2<f64>,
        B: Array2<f64>,
        C: Array2<f64>,
        f: Array2<f64>,
        v: Array2<f64>,
        W3: Array2<f64>,
        W4: Array2<f64>,
        desired_control_traj: Array2<f64>,
        current_timestep: i64,
    }

    impl Controller {
        fn form_lifted_matrices(&self) {}

        fn propagate_dynamics(&self, control_input: &Array2<f64>, state: &Array2<f64>) {}

        fn compute_control_inputs(&self) {}
    }
}

use simulator::system_simulate;
use std::env;

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
    let (y_test, x_test) = system_simulate(&A, &B, &C, &input_test, &x0_test);

    // Draw the response
    let root = BitMapBackend::new("step_response.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_y = y_test.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y = y_test.iter().cloned().fold(f64::INFINITY, f64::min);
    let mut chart = ChartBuilder::on(&root)
        .caption("System Output Y", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0..y_test.shape()[1] as i32, min_y..max_y)?;

    chart.configure_mesh().draw()?;

    // Plot input
    let series_input: Vec<(i32, f64)> = input_test
        .row(0)
        .iter()
        .enumerate()
        .map(|(i, &val)| (i as i32, val as f64))
        .collect();

    chart
        .draw_series(LineSeries::new(series_input, &Palette99::pick(0)))?
        .label(format!("Output {}", 0))
        .legend(move |(x, y)| PathElement::new([(x, y), (x + 20, y)], &Palette99::pick(0)));

    // Plot system response
    let series_y: Vec<(i32, f64)> = y_test
        .row(0)
        .iter()
        .enumerate()
        .map(|(i, &val)| (i as i32, val as f64))
        .collect();

    chart.draw_series(LineSeries::new(series_y, &Palette99::pick(1)))?;

    chart
        .configure_series_labels()
        .background_style(&WHITE)
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}
