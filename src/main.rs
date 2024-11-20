mod simulator {
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
                mat_y.slice_mut(s![.., i]).assign(&(mat_x.dot(x0)));
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
}

mod controller {
    use ndarray::{array, s, Array1, Array2, Axis};
    use ndarray_linalg::Inverse;

    pub struct Controller {
        mat_a: Array2<f64>,
        mat_b: Array2<f64>,
        mat_c: Array2<f64>,
        f: usize,
        v: usize,
        mat_w3: Array2<f64>,
        mat_w4: Array2<f64>,
        desired_ctrl_traj_total: Array2<f64>,
        current_timestep: usize,
        mat_o: Array2<f64>,
        mat_m: Array2<f64>,
        gain_matrix: Array2<f64>,
        states: Array2<f64>,
        outputs: Array2<f64>,
        inputs: Array2<f64>,
    }

    impl Controller {
        pub fn form_lifted_matrices(
            mat_a: &Array2<f64>,
            mat_b: &Array2<f64>,
            mat_c: &Array2<f64>,
            f: usize,
            v: usize,
            mat_w3: &Array2<f64>,
            mat_w4: &Array2<f64>,
        ) -> Result<(Array2<f64>, Array2<f64>, Array2<f64>), Box<dyn std::error::Error>> {
            let n = mat_a.nrows();
            let r = mat_c.nrows();
            let m = mat_b.ncols();

            // Lifted matrix O
            let mut mat_o: Array2<f64> = Array2::zeros((f * r, n));

            let mut pow_a = mat_a.clone();
            for i in 0..f {
                if i != 0 {
                    pow_a.assign(&(pow_a.dot(mat_a)));
                }

                mat_o
                    .slice_mut(s![i * r..(i + 1) * r, ..])
                    .assign(&(mat_c.dot(&pow_a)));
            }

            // Lifted matrix M
            let mut mat_m: Array2<f64> = Array2::zeros((f * r, v * m));

            for i in 0..f {
                // Until the control horizon
                if i < v {
                    for j in 0..(i + 1) {
                        if j == 0 {
                            pow_a = Array2::eye(n);
                        } else {
                            pow_a.assign(&(pow_a.dot(mat_a)));
                        }

                        mat_m
                            .slice_mut(s![i * r..(i + 1) * r, (i - j) * m..(i - j + 1) * m])
                            .assign(&(mat_c.dot(&pow_a).dot(mat_b)));
                    }
                } else {
                    for j in 0..v {
                        // Here we form the last entry
                        if j == 0 {
                            let mut sum_last: Array2<f64> = Array2::zeros((n, n));
                            for s in 0..i - v + 2 {
                                if s == 0 {
                                    pow_a = Array2::eye(n);
                                } else {
                                    pow_a.assign(&(pow_a.dot(mat_a)));
                                }

                                sum_last = sum_last + &pow_a;
                            }

                            mat_m
                                .slice_mut(s![i * r..(i + 1) * r, (v - 1) * m..(v) * m])
                                .assign(&(mat_c.dot(&sum_last).dot(mat_b)));
                        } else {
                            pow_a.assign(&(pow_a.dot(mat_a)));

                            mat_m
                                .slice_mut(s![i * r..(i + 1) * r, (v - 1 - j) * m..(v - j) * m])
                                .assign(&(mat_c.dot(&pow_a).dot(mat_b)));
                        }
                    }
                }
            }

            let tmp1 = mat_m.t().dot(mat_w4).dot(&mat_m);
            let tmp2: Array2<f64> = (tmp1 + mat_w3).to_owned().inv()?;
            let gain_matrix = tmp2.dot(&mat_m.t()).dot(mat_w4);

            Ok((mat_o, mat_m, gain_matrix))
        }

        fn propagate_dynamics(
            &self,
            control_input: &Array2<f64>,
            state: &Array1<f64>,
        ) -> (Array2<f64>, Array2<f64>) {
            let mut x_kp1 = Array2::zeros((self.mat_a.shape()[0], 1));
            let mut y_k = Array2::zeros((self.mat_c.shape()[0], 1));

            x_kp1.assign(&(self.mat_a.dot(state) + self.mat_b.dot(control_input)));
            y_k.assign(&(self.mat_c.dot(state)));

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
                    .mat_o
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
            let _ = self.states.append(Axis(0), (&state_kp1).into());
            let _ = self.outputs.append(Axis(0), (&output_k).into());
            let _ = self.inputs.append(Axis(0), (&input_applied).into());

            // Increment the time step
            self.current_timestep = self.current_timestep + 1;
        }

        fn new(
            mat_a: &Array2<f64>,
            mat_b: &Array2<f64>,
            mat_c: &Array2<f64>,
            f: usize,
            v: usize,
            mat_w3: &Array2<f64>,
            mat_w4: &Array2<f64>,
            x0: Array1<f64>,
            desired_ctrl_traj: &Array2<f64>,
        ) -> Result<Controller, Box<dyn std::error::Error>> {
            // Form the lifted system matrices and vectors
            // the gain matrix is used to compute the solution
            // here we pre-compute it to save computational time
            let (mat_o, mat_m, gain_matrix) =
                Self::form_lifted_matrices(mat_a, mat_b, mat_c, f, v, mat_w3, mat_w4)?;

            // We store the state vectors of the controlled state trajectory
            let mut states: Array2<f64> = Array2::zeros((1, 2));
            states.slice_mut(s![1, ..]).assign(&x0);

            // We store the computed inputs
            let inputs: Array2<f64> = array![[]];

            // # we store the output vectors of the controlled state trajectory
            let mut outputs: Array2<f64> = Array2::zeros((1, 2));
            outputs.slice_mut(s![1, ..]).assign(&(mat_c.dot(&x0)));

            Ok(Controller {
                mat_a: mat_a.clone(),
                mat_b: mat_b.clone(),
                mat_c: mat_c.clone(),
                f: f,
                v: v,
                mat_w3: mat_w3.clone(),
                mat_w4: mat_w4.clone(),
                desired_ctrl_traj_total: desired_ctrl_traj.clone(),
                current_timestep: 0,
                mat_o: mat_o,
                mat_m: mat_m,
                gain_matrix: gain_matrix,
                states: states,
                inputs: inputs,
                outputs: outputs,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_form_lifted_matrix() -> Result<(), Box<dyn std::error::Error>> {
        // We want to do a simple test with a state space of two and one control variable. We first start with a one step horizon
        let A: Array2<f64> = array![[1.0, 0.0], [0.0, 2.0]];
        let B: Array2<f64> = array![[3.0], [4.0]];
        let C: Array2<f64> = array![[5.0, 6.0]];
        let f = 3usize; // Prediction horizon
        let v = 3usize; // Control horizon
        let W3: Array2<f64> = array![[5.0, -3.0, 0.0], [-3.0, 7.0, -4.0], [0.0, -4.0, 4.0]];
        let W4: Array2<f64> = array![[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [7.0, 0.0, 0.0]];

        let (O, M, gain_matrix) =
            controller::Controller::form_lifted_matrices(&A, &B, &C, f, v, &W3, &W4)?;

        // For a simple test case with horizon of 1, O is equal to A
        // assert_eq!(O.shape(), &[f * C.nrows(), A.nrows()]);
        // assert_eq!(O, A);

        // For a simple test case with horizon of 1, O is equal to A
        let expected_M = array![[39.0, 0.0, 0.0], [63.0, 39.0, 0.0], [111.0, 63.0, 39.0]];
        assert_eq!(M, expected_M);

        Ok(())
    }
}

use ndarray::{array, s, Array, Array1, Array2};
use ndarray_linalg::{Eig, Inverse};
use plotters::prelude::*;
use simulator::system_simulate;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env::set_var("RUST_BACKTRACE", "1");

    // Define parameters
    let f = 20usize;
    let v = 20usize;
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

    let r = 1usize;
    let m = 1usize; // number of inputs and outputs
    let _n = 4usize; // state dimension

    // Discretization constant
    let sampling = 0.05;

    // Model discretization
    let mat_i: Array2<f64> = Array::eye(mat_ac.nrows());
    let mut mat_a: Array2<f64> = mat_i - sampling * mat_ac.clone();
    mat_a = mat_a.inv()?;
    let mat_b = mat_a.dot(&(sampling * mat_bc));
    let mat_c = mat_cc;

    // check the eigenvalues
    let _eigen_a = mat_ac.eig()?;
    let _eigen_aid = mat_a.eig()?;

    let time_sample_test = 200;

    // Compute the system's step response
    let input_test = 10.0 * Array2::ones((1, time_sample_test));
    let x0_test = Array1::zeros(4);

    // // # simulate the discrete-time system
    let (y_test, _x_test) = system_simulate(&mat_a, &mat_b, &mat_c, &input_test, &x0_test);

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
        .build_cartesian_2d(0..y_test.ncols() as i32, min_y..max_y)?;

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

    let mut mat_w1: Array2<f64> = Array2::zeros((v * m, v * m));

    for i in 0..v {
        if i == 0 {
            mat_w1
                .slice_mut(s![i * m..(i + 1) * m, i * m..(i + 1) * m])
                .assign(&Array2::eye(m));
        } else {
            mat_w1
                .slice_mut(s![i * m..(i + 1) * m, i * m..(i + 1) * m])
                .assign(&Array2::eye(m));
            mat_w1
                .slice_mut(s![i * m..(i + 1) * m, (i - 1) * m..(i) * m])
                .assign(&Array2::eye(m));
        }
    }

    // mat_w2 matrix
    let mat_q0 = array![0.0000000011f64];
    let math_q_other = array![0.0001f64];

    let mut mat_w2: Array2<f64> = Array2::zeros((v * m, v * m));

    for i in 0..v {
        if i == 0 {
            mat_w2
                .slice_mut(s![i * m..(i + 1) * m, i * m..(i + 1) * m])
                .assign(&mat_q0);
        } else {
            mat_w2
                .slice_mut(s![i * m..(i + 1) * m, i * m..(i + 1) * m])
                .assign(&math_q_other);
        }
    }

    // W3 matrix
    let _mat_w3 = mat_w1.t().dot(&(mat_w2.dot(&mat_w1)));

    // W4 matrix
    let mut mat_w4: Array2<f64> = Array2::zeros((f * r, f * r));

    let pred_weight = array![10f64];

    for i in 0..f {
        mat_w4
            .slice_mut(s![i * r..(i + 1) * r, i * r..(i + 1) * r])
            .assign(&pred_weight);
    }

    let time_steps = 300;

    // Define a step trajectory
    let desired_traj : Array2<f64> = 0.3 * Array2::ones((time_steps, 1));
 
    Ok(())
}
