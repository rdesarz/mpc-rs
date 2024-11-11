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
    use ndarray::{s, Array1, Array2, Axis};
    use ndarray_linalg::Inverse;

    pub struct Controller {
        A: Array2<f64>,
        B: Array2<f64>,
        C: Array2<f64>,
        f: usize,
        v: usize,
        W3: Array2<f64>,
        W4: Array2<f64>,
        desired_ctrl_traj_total: Array2<f64>,
        current_timestep: usize,
        O: Array2<f64>,
        M: Array2<f64>,
        gain_matrix: Array2<f64>,
        states: Array2<f64>,
        outputs: Array2<f64>,
        inputs: Array2<f64>,
    }

    impl Controller {
        // fn new(A,B,C,f,v,W3,W4,x0,desiredControlTrajectoryTotal) -> Controller {
        //     Point { x: x, y: y }
        // }

        fn form_lifted_matrices(
            &self,
        ) -> Result<(Array2<f64>, Array2<f64>, Array2<f64>), Box<dyn std::error::Error>> {
            let n = self.A.nrows();
            let r = self.C.nrows();
            let m = self.B.ncols();

            // Lifted matrix O
            let mut O: Array2<f64> = Array2::zeros((self.f * r, n));

            let mut powA = self.A.clone();
            for i in 0..self.f {
                if i != 0 {
                    powA.assign(&(powA.dot(&self.A)));
                }

                O.slice_mut(s![i * r..(i + 1) * r, ..])
                    .assign(&(self.C.dot(&powA)));
            }

            // Lifted matrix M
            let mut M: Array2<f64> = Array2::zeros((self.f * r, self.v * m));

            for i in 0..self.f {
                // Until the control horizon
                if i < self.v {
                    for j in 0..(i + 1) {
                        if j == 0 {
                            powA = Array2::eye(n);
                        } else {
                            powA.assign(&(powA.dot(&self.A)));
                        }

                        M.slice_mut(s![i * r..(i + 1) * r, (i - j) * m..(i - j + 1) * m])
                            .assign(&(self.C.dot(&powA).dot(&self.B)));
                    }
                } else {
                    for j in 0..self.v {
                        // Here we form the last entry
                        if j == 0 {
                            let mut sum_last: Array2<f64> = Array2::zeros((n, n));
                            for s in 0..i - self.v + 2 {
                                if s == 0 {
                                    powA = Array2::eye(n);
                                } else {
                                    powA.assign(&(powA.dot(&self.A)));
                                }

                                sum_last = sum_last + &powA;
                            }

                            M.slice_mut(s![i * r..(i + 1) * r, (self.v - 1) * m..(self.v) * m])
                                .assign(&(self.C.dot(&sum_last).dot(&self.B)));
                        } else {
                            powA.assign(&(powA.dot(&self.A)));

                            M.slice_mut(s![
                                i * r..(i + 1) * r,
                                (self.v - 1 - j) * m..(self.v - j) * m
                            ])
                            .assign(&(self.C.dot(&powA).dot(&self.B)));
                        }
                    }
                }
            }

            let tmp1 = M.t().dot(&self.W4).dot(&M);
            let tmp2: Array2<f64> = (tmp1 + &self.W3).to_owned().inv()?;
            let gain_matrix = tmp2.dot(&M.t()).dot(&self.W4);

            Ok((O, M, gain_matrix))
        }

        fn propagate_dynamics(
            &self,
            control_input: &Array2<f64>,
            state: &Array1<f64>,
        ) -> (Array2<f64>, Array2<f64>) {
            let mut x_kp1 = Array2::zeros((self.A.shape()[0], 1));
            let mut y_k = Array2::zeros((self.C.shape()[0], 1));

            x_kp1.assign(&(self.A.dot(state) + self.B.dot(control_input)));
            y_k.assign(&(self.C.dot(state)));

            (x_kp1, y_k)
        }

        fn compute_control_inputs(&mut self) {
            // Extract the segment of the desired control trajectory
            let desired_ctrl_traj = self
                .desired_ctrl_traj_total
                .slice(s![
                    self.current_timestep..(self.current_timestep + self.f),
                    ..
                ])
                .to_owned();

            // Compute the vector s
            let vec_s = desired_ctrl_traj
                - self
                    .O
                    .dot(&self.states.slice(s![self.current_timestep, ..]));

            // Compute the control sequence
            let input_sequence_computed = self.gain_matrix.dot(&vec_s);
            let mut input_applied: Array2<f64> = Array2::zeros((1, 1));
            input_applied[[0, 0]] = input_sequence_computed[[0, 0]];

            // Compute the next state and output
            let (state_kp1, output_k) = self.propagate_dynamics(
                &input_applied,
                &self.states.slice(s![self.current_timestep, ..]).to_owned(),
            );

            // Append the lists
            self.states.append(Axis(0), (&state_kp1).into());
            self.outputs.append(Axis(0), (&output_k).into());
            self.inputs.append(Axis(0), (&input_applied).into());

            // Increment the time step
            self.current_timestep = self.current_timestep + 1;
        }
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
