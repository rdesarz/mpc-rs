
use nalgebra as na;

pub struct Controller {
    mat_a: na::DMatrix<f64>,
    mat_b: na::DMatrix<f64>,
    mat_c: na::DMatrix<f64>,
    f: usize,
    current_timestep: usize,
    mat_o: na::DMatrix<f64>,
    gain_matrix: na::DMatrix<f64>,
    pub states: na::DMatrix<f64>,
    pub desired_ctrl_traj_total: na::DMatrix<f64>,
    pub outputs: na::DMatrix<f64>,
    pub inputs: na::DMatrix<f64>,
}

impl Controller {
    // pub fn form_lifted_matrices(
    //     mat_a: &na::DMatrix<f64>,
    //     mat_b: &na::DMatrix<f64>,
    //     mat_c: &na::DMatrix<f64>,
    //     f: usize,
    //     v: usize,
    //     mat_w3: &na::DMatrix<f64>,
    //     mat_w4: &na::DMatrix<f64>,
    // ) -> Result<(na::DMatrix<f64>, na::DMatrix<f64>, na::DMatrix<f64>), Box<dyn std::error::Error>> {
    //     let n = mat_a.nrows();
    //     let r = mat_c.nrows();
    //     let m = mat_b.ncols();

    //     // Lifted matrix O
    //     let mut mat_o: na::DMatrix<f64> = na::DMatrix<f64>::zeros(f * r, n);

    //     let mut pow_a = mat_a.clone();
    //     for i in 0..f {
    //         if i != 0 {
    //             pow_a = pow_a * &mat_a;
    //         }

    //         mat_o
    //             .slice_mut(s![i * r..(i + 1) * r, ..])
    //             .assign(&(mat_c.dot(&pow_a)));
    //     }

    //     // Lifted matrix M
    //     let mut mat_m: na::DMatrix<f64> = Array2::zeros((f * r, v * m));

    //     for i in 0..f {
    //         // Until the control horizon
    //         if i < v {
    //             for j in 0..(i + 1) {
    //                 if j == 0 {
    //                     pow_a = Array2::eye(n);
    //                 } else {
    //                     pow_a.assign(&(pow_a.dot(mat_a)));
    //                 }

    //                 mat_m
    //                     .slice_mut(s![i * r..(i + 1) * r, (i - j) * m..(i - j + 1) * m])
    //                     .assign(&(mat_c.dot(&(pow_a.dot(mat_b)))));
    //             }
    //         } else {
    //             for j in 0..v {
    //                 // Here we form the last entry
    //                 if j == 0 {
    //                     let mut sum_last: na::DMatrix<f64> = Array2::zeros((n, n));
    //                     for s in 0..i - v + 2 {
    //                         if s == 0 {
    //                             pow_a = Array2::eye(n);
    //                         } else {
    //                             pow_a.assign(&(pow_a.dot(mat_a)));
    //                         }

    //                         sum_last = sum_last + &pow_a;
    //                     }

    //                     mat_m
    //                         .slice_mut(s![i * r..(i + 1) * r, (v - 1) * m..(v) * m])
    //                         .assign(&(mat_c.dot(&(sum_last.dot(mat_b)))));
    //                 } else {
    //                     pow_a.assign(&(pow_a.dot(mat_a)));

    //                     mat_m
    //                         .slice_mut(s![i * r..(i + 1) * r, (v - 1 - j) * m..(v - j) * m])
    //                         .assign(&(mat_c.dot(&(pow_a.dot(mat_b)))));
    //                 }
    //             }
    //         }
    //     }

    //     let tmp1 = mat_m.t().dot(&(mat_w4.dot(&mat_m)));
    //     let tmp2: na::DMatrix<f64> = (tmp1 + mat_w3).to_owned().inv()?;
    //     let gain_matrix = tmp2.dot(&(mat_m.t().dot(mat_w4)));

    //     Ok((mat_o, mat_m, gain_matrix))
    // }

    pub fn propagate_dynamics(
        &self,
        control_input: &na::DVector<f64>,
        state: &na::DVector<f64>,
    ) -> (na::DVector<f64>, na::DVector<f64>) {
        let mut x_kp1 = na::DVector::zeros(self.mat_a.nrows());
        let mut y_k = na::DVector::zeros(self.mat_c.nrows());

        x_kp1.copy_from(&(&self.mat_a * state + &self.mat_b * control_input));
        y_k.copy_from(&(&self.mat_c * state));

        (x_kp1, y_k)
    }

    // pub fn compute_control_inputs(&mut self) {
    //     // Extract the segment of the desired control trajectory
    //     let desired_ctrl_traj = self
    //         .desired_ctrl_traj_total
    //         .slice(s![
    //             self.current_timestep..(self.current_timestep + self.f),
    //             ..
    //         ])
    //         .to_owned();

    //     // Compute the vector s
    //     let vec_s = (desired_ctrl_traj.t().to_owned()
    //         - self
    //             .mat_o
    //             .dot(&self.states.slice(s![self.current_timestep, ..])))
    //     .t()
    //     .to_owned();

    //     // Compute the control sequence
    //     let input_sequence_computed = self.gain_matrix.dot(&vec_s);
    //     let mut input_applied: Array1<f64> = Array1::zeros(1);
    //     input_applied[0] = input_sequence_computed[[0, 0]];

    //     // Compute the next state and output
    //     let (state_kp1, output_k) = self.propagate_dynamics(
    //         &input_applied,
    //         &self.states.slice(s![self.current_timestep, ..]).to_owned(),
    //     );

    //     // Append the lists
    //     self.states = ndarray::concatenate(
    //         Axis(0),
    //         &[self.states.view(), state_kp1.insert_axis(Axis(0)).view()],
    //     )
    //     .unwrap();

    //     self.outputs = ndarray::concatenate(
    //         Axis(0),
    //         &[self.outputs.view(), output_k.insert_axis(Axis(0)).view()],
    //     )
    //     .unwrap();

    //     self.inputs = ndarray::concatenate(
    //         Axis(1),
    //         &[
    //             self.inputs.view(),
    //             input_applied.insert_axis(Axis(0)).view(),
    //         ],
    //     )
    //     .unwrap();

    //     // Increment the time step
    //     self.current_timestep = self.current_timestep + 1;
    // }

    // pub fn new(
    //     mat_a: na::DMatrix<f64>,
    //     mat_b: na::DMatrix<f64>,
    //     mat_c: na::DMatrix<f64>,
    //     f: usize,
    //     v: usize,
    //     mat_w3: &na::DMatrix<f64>,
    //     mat_w4: &na::DMatrix<f64>,
    //     x0: Array1<f64>,
    //     desired_ctrl_traj: &na::DMatrix<f64>,
    // ) -> Result<Controller, Box<dyn std::error::Error>> {
    //     // Form the lifted system matrices and vectors
    //     // the gain matrix is used to compute the solution
    //     // here we pre-compute it to save computational time
    //     let (mat_o, _, gain_matrix) =
    //         Self::form_lifted_matrices(&mat_a, &mat_b, &mat_c, f, v, mat_w3, mat_w4)?;

    //     // We store the state vectors of the controlled state trajectory
    //     let mut states: na::DMatrix<f64> = Array2::zeros((1, x0.shape()[0]));
    //     states.slice_mut(s![0, ..]).assign(&x0);

    //     // We store the computed inputs
    //     let inputs: na::DMatrix<f64> = array![[]];

    //     // # we store the output vectors of the controlled state trajectory
    //     let mut outputs: na::DMatrix<f64> = Array2::zeros((1, mat_c.nrows()));
    //     outputs.slice_mut(s![0, ..]).assign(&(mat_c.dot(&x0)));

    //     Ok(Controller {
    //         mat_a: mat_a,
    //         mat_b: mat_b,
    //         mat_c: mat_c,
    //         f: f,
    //         desired_ctrl_traj_total: desired_ctrl_traj.clone(),
    //         current_timestep: 0,
    //         mat_o: mat_o,
    //         gain_matrix: gain_matrix,
    //         states: states,
    //         inputs: inputs,
    //         outputs: outputs,
    //     })
    // }
}
