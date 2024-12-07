
use ndarray::{array, s, Array1, Array2, Axis};
use ndarray_linalg::Inverse;

pub struct Controller {
    mat_a: Array2<f64>,
    mat_b: Array2<f64>,
    mat_c: Array2<f64>,
    f: usize,
    current_timestep: usize,
    mat_o: Array2<f64>,
    gain_matrix: Array2<f64>,
    pub states: Array2<f64>,
    pub desired_ctrl_traj_total: Array2<f64>,
    pub outputs: Array2<f64>,
    pub inputs: Array2<f64>,
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
                        .assign(&(mat_c.dot(&(pow_a.dot(mat_b)))));
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
                            .assign(&(mat_c.dot(&(sum_last.dot(mat_b)))));
                    } else {
                        pow_a.assign(&(pow_a.dot(mat_a)));

                        mat_m
                            .slice_mut(s![i * r..(i + 1) * r, (v - 1 - j) * m..(v - j) * m])
                            .assign(&(mat_c.dot(&(pow_a.dot(mat_b)))));
                    }
                }
            }
        }

        let tmp1 = mat_m.t().dot(&(mat_w4.dot(&mat_m)));
        let tmp2: Array2<f64> = (tmp1 + mat_w3).to_owned().inv()?;
        let gain_matrix = tmp2.dot(&(mat_m.t().dot(mat_w4)));

        Ok((mat_o, mat_m, gain_matrix))
    }

    pub fn propagate_dynamics(
        &self,
        control_input: &Array1<f64>,
        state: &Array1<f64>,
    ) -> (Array1<f64>, Array1<f64>) {
        let mut x_kp1 = Array1::zeros(self.mat_a.nrows());
        let mut y_k = Array1::zeros(self.mat_c.nrows());

        x_kp1.assign(&(self.mat_a.dot(state) + self.mat_b.dot(control_input)));
        y_k.assign(&(self.mat_c.dot(state)));

        (x_kp1, y_k)
    }

    pub fn compute_control_inputs(&mut self) {
        // Extract the segment of the desired control trajectory
        let desired_ctrl_traj = self
            .desired_ctrl_traj_total
            .slice(s![
                self.current_timestep..(self.current_timestep + self.f),
                ..
            ])
            .to_owned();

        // Compute the vector s
        let vec_s = (desired_ctrl_traj.t().to_owned()
            - self
                .mat_o
                .dot(&self.states.slice(s![self.current_timestep, ..])))
        .t()
        .to_owned();

        // Compute the control sequence
        let input_sequence_computed = self.gain_matrix.dot(&vec_s);
        let mut input_applied: Array1<f64> = Array1::zeros(1);
        input_applied[0] = input_sequence_computed[[0, 0]];

        // Compute the next state and output
        let (state_kp1, output_k) = self.propagate_dynamics(
            &input_applied,
            &self.states.slice(s![self.current_timestep, ..]).to_owned(),
        );

        // Append the lists
        self.states = ndarray::concatenate(
            Axis(0),
            &[self.states.view(), state_kp1.insert_axis(Axis(0)).view()],
        )
        .unwrap();

        self.outputs = ndarray::concatenate(
            Axis(0),
            &[self.outputs.view(), output_k.insert_axis(Axis(0)).view()],
        )
        .unwrap();

        self.inputs = ndarray::concatenate(
            Axis(1),
            &[
                self.inputs.view(),
                input_applied.insert_axis(Axis(0)).view(),
            ],
        )
        .unwrap();

        // Increment the time step
        self.current_timestep = self.current_timestep + 1;
    }

    pub fn new(
        mat_a: Array2<f64>,
        mat_b: Array2<f64>,
        mat_c: Array2<f64>,
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
        let (mat_o, _, gain_matrix) =
            Self::form_lifted_matrices(&mat_a, &mat_b, &mat_c, f, v, mat_w3, mat_w4)?;

        // We store the state vectors of the controlled state trajectory
        let mut states: Array2<f64> = Array2::zeros((1, x0.shape()[0]));
        states.slice_mut(s![0, ..]).assign(&x0);

        // We store the computed inputs
        let inputs: Array2<f64> = array![[]];

        // # we store the output vectors of the controlled state trajectory
        let mut outputs: Array2<f64> = Array2::zeros((1, mat_c.nrows()));
        outputs.slice_mut(s![0, ..]).assign(&(mat_c.dot(&x0)));

        Ok(Controller {
            mat_a: mat_a,
            mat_b: mat_b,
            mat_c: mat_c,
            f: f,
            desired_ctrl_traj_total: desired_ctrl_traj.clone(),
            current_timestep: 0,
            mat_o: mat_o,
            gain_matrix: gain_matrix,
            states: states,
            inputs: inputs,
            outputs: outputs,
        })
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_form_lifted_matrix() -> Result<(), Box<dyn std::error::Error>> {
        // We want to do a simple test with a state space of two and one control variable. We first start with a one step horizon
        let mat_a: Array2<f64> = array![[1.0, 0.0], [0.0, 2.0]];
        let mat_b: Array2<f64> = array![[3.0], [4.0]];
        let mat_c: Array2<f64> = array![[5.0, 6.0]];
        let f = 3usize; // Prediction horizon
        let v = 3usize; // Control horizon
        let mat_w3: Array2<f64> = array![[5.0, -3.0, 0.0], [-3.0, 7.0, -4.0], [0.0, -4.0, 4.0]];
        let mat_w4: Array2<f64> = array![[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [7.0, 0.0, 0.0]];

        let (mat_o, mat_m, gain_matrix) = controller::Controller::form_lifted_matrices(
            &mat_a, &mat_b, &mat_c, f, v, &mat_w3, &mat_w4,
        )?;

        let expected_mat_o = array![[5.0, 12.0], [5.0, 24.0], [5.0, 48.0]];
        assert_eq!(mat_o, expected_mat_o);

        let expected_mat_m = array![[39.0, 0.0, 0.0], [63.0, 39.0, 0.0], [111.0, 63.0, 39.0]];
        assert_eq!(mat_m, expected_mat_m);

        let expected_gain_matrix = array![
            [0.02564045344996876, 0.0, 0.0],
            [0.03269213603500207, 0.0, 0.0],
            [0.034215165580676776, 0.0, 0.0]
        ];
        assert_eq!(gain_matrix, expected_gain_matrix);

        Ok(())
    }

    #[test]
    fn test_controller_creation() -> Result<(), Box<dyn std::error::Error>> {
        // Define a step trajectory
        let desired_traj: Array2<f64> = 0.3 * Array2::ones((100, 1));

        // Set the initial state
        let x0 = Array::ones(2);

        // We want to do a simple test with a state space of two and one control variable. We first start with a one step horizon
        let mat_a: Array2<f64> = array![[1.0, 0.0], [0.0, 2.0]];
        let mat_b: Array2<f64> = array![[3.0], [4.0]];
        let mat_c: Array2<f64> = array![[5.0, 6.0]];
        let f = 3usize; // Prediction horizon
        let v = 3usize; // Control horizon
        let mat_w3: Array2<f64> = array![[5.0, -3.0, 0.0], [-3.0, 7.0, -4.0], [0.0, -4.0, 4.0]];
        let mat_w4: Array2<f64> = array![[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [7.0, 0.0, 0.0]];

        // Create the controller
        let mpc = Controller::new(
            mat_a,
            mat_b,
            mat_c,
            f,
            v,
            &mat_w3,
            &mat_w4,
            x0,
            &desired_traj,
        )?;

        // State is equal to initial state
        assert_eq!(mpc.states, array![[1.0, 1.0]]);

        // No input stored after initialisation
        assert_eq!(mpc.inputs, array![[]]);

        // First output is state times output matrix
        assert_eq!(mpc.outputs, array![[11.0]]);

        Ok(())
    }

    #[test]
    fn test_compute_control_inputs() -> Result<(), Box<dyn std::error::Error>> {
        // Define a step trajectory
        let desired_traj: Array2<f64> = 0.3 * Array2::ones((100, 1));

        // Set the initial state
        let x0 = Array::ones(2);

        // We want to do a simple test with a state space of two and one control variable. We first start with a one step horizon
        let mat_a: Array2<f64> = array![[1.0, 0.0], [0.0, 2.0]];
        let mat_b: Array2<f64> = array![[3.0], [4.0]];
        let mat_c: Array2<f64> = array![[5.0, 6.0]];
        let f = 3usize; // Prediction horizon
        let v = 3usize; // Control horizon
        let mat_w3: Array2<f64> = array![[5.0, -3.0, 0.0], [-3.0, 7.0, -4.0], [0.0, -4.0, 4.0]];
        let mat_w4: Array2<f64> = array![[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [7.0, 0.0, 0.0]];

        // Create the controller
        let mut mpc = Controller::new(
            mat_a,
            mat_b,
            mat_c,
            f,
            v,
            &mat_w3,
            &mat_w4,
            x0,
            &desired_traj,
        )?;

        mpc.compute_control_inputs();

        assert_eq!(mpc.outputs, array![[11.0], [11.0]]);

        Ok(())
    }
}
