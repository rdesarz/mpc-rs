use nalgebra as na;

use crate::mpc::linear_discrete_model::LinearDiscreteModel;

pub struct Controller<'a> {
    model: &'a dyn LinearDiscreteModel,
    f: usize,
    current_timestep: usize,
    mat_o: na::DMatrix<f64>,
    gain_matrix: na::DMatrix<f64>,
    pub states: na::DMatrix<f64>,
    pub desired_ctrl_traj_total: na::DMatrix<f64>,
    pub outputs: na::DMatrix<f64>,
    pub inputs: na::DMatrix<f64>,
}

impl<'a> Controller<'a> {
    pub fn form_lifted_matrices(
        mat_a: &na::DMatrix<f64>,
        mat_b: &na::DMatrix<f64>,
        mat_c: &na::DMatrix<f64>,
        f: usize,
        v: usize,
        mat_w3: &na::DMatrix<f64>,
        mat_w4: &na::DMatrix<f64>,
    ) -> Result<(na::DMatrix<f64>, na::DMatrix<f64>, na::DMatrix<f64>), Box<dyn std::error::Error>>
    {
        let n = mat_a.nrows();
        let r = mat_c.nrows();
        let m = mat_b.ncols();

        // Lifted matrix O
        let mut mat_o = na::DMatrix::<f64>::zeros(f * r, n);

        let mut pow_a = mat_a.clone();
        for i in 0..f {
            if i != 0 {
                pow_a = &pow_a * mat_a;
            }

            mat_o
                .view_mut((i * r, 0), (r, n))
                .copy_from(&(&*mat_c * &pow_a));
        }

        // Lifted matrix M
        let mut mat_m: na::DMatrix<f64> = na::DMatrix::zeros(f * r, v * m);

        for i in 0..f {
            // Until the control horizon
            if i < v {
                for j in 0..=i {
                    if j == 0 {
                        pow_a = na::DMatrix::identity(n, n);
                    } else {
                        pow_a = &pow_a * mat_a;
                    }

                    mat_m
                        .view_mut((i * r, (i - j) * m), (r, m))
                        .copy_from(&(&*mat_c * &pow_a * mat_b));
                }
            } else {
                for j in 0..v {
                    // Here we form the last entry
                    if j == 0 {
                        let mut sum_last: na::DMatrix<f64> = na::DMatrix::zeros(n, n);
                        for s in 0..(i - v + 2) {
                            if s == 0 {
                                pow_a = na::DMatrix::identity(n, n);
                            } else {
                                pow_a = &pow_a * mat_a;
                            }

                            sum_last += &pow_a;
                        }

                        mat_m
                            .view_mut((i * r, (v - 1) * m), (r, m))
                            .copy_from(&(&*mat_c * &sum_last * mat_b));
                    } else {
                        pow_a = &pow_a * mat_a;

                        mat_m
                            .view_mut((i * r, (v - 1 - j) * m), (r, m))
                            .copy_from(&(&*mat_c * &pow_a * mat_b));
                    }
                }
            }
        }

        let tmp1 = mat_m.transpose() * mat_w4 * &mat_m;
        let tmp2: na::DMatrix<f64> = (tmp1 + mat_w3)
            .try_inverse()
            .ok_or("Matrix inversion failed")?;
        let gain_matrix = tmp2 * mat_m.transpose() * mat_w4;

        Ok((mat_o, mat_m, gain_matrix))
    }

    pub fn propagate_dynamics(
        &self,
        control_input: &na::DVector<f64>,
        state: &na::DVector<f64>,
    ) -> (na::DVector<f64>, na::DVector<f64>) {
        let mut x_kp1 = na::DVector::zeros(self.model.get_mat_a().nrows());
        let mut y_k = na::DVector::zeros(self.model.get_mat_c().nrows());

        x_kp1.copy_from(&(self.model.get_mat_a() * state + self.model.get_mat_b() * control_input));
        y_k.copy_from(&(self.model.get_mat_c() * state));

        (x_kp1, y_k)
    }

    pub fn compute_control_inputs(&mut self) {
        // Extract the segment of the desired control trajectory
        let desired_ctrl_traj = self
            .desired_ctrl_traj_total
            .view_range(self.current_timestep..self.current_timestep + self.f, ..)
            .into_owned();

        // Compute the vector s
        let vec_s = (desired_ctrl_traj - &self.mat_o * self.states.column(self.current_timestep))
            .into_owned();

        // Compute the control sequence
        let input_sequence_computed = &self.gain_matrix * vec_s;
        let mut input_applied = na::DVector::<f64>::zeros(1);
        input_applied[0] = input_sequence_computed[(0, 0)];

        // Compute the next state and output
        let (state_kp1, output_k) = self.propagate_dynamics(
            &input_applied,
            &self.states.column(self.current_timestep).into_owned(),
        );

        // Append the lists
        self.states = na::stack![self.states, state_kp1];
        self.outputs = na::stack![self.outputs, output_k];

        if self.inputs.shape() == (0, 0) {
            self.inputs.resize_mut(1, input_applied.nrows(), 0.);
            self.inputs.view_range_mut(0, ..).copy_from(&input_applied);
        } else {
            self.inputs = na::stack![self.inputs, input_applied];
        }

        // Increment the time step
        self.current_timestep = self.current_timestep + 1;
    }

    pub fn new(
        model: &'a dyn LinearDiscreteModel,
        f: usize,
        v: usize,
        mat_w3: &na::DMatrix<f64>,
        mat_w4: &na::DMatrix<f64>,
        x0: na::DVector<f64>,
        desired_ctrl_traj: &na::DMatrix<f64>,
    ) -> Result<Controller<'a>, Box<dyn std::error::Error>> {
        // Form the lifted system matrices and vectors
        // the gain matrix is used to compute the solution
        // here we pre-compute it to save computational time
        let (mat_o, _, gain_matrix) =
            Self::form_lifted_matrices(model.get_mat_a(), model.get_mat_b(), model.get_mat_c(), f, v, mat_w3, mat_w4)?;

        // We store the state vectors of the controlled state trajectory. States are stored as column
        let mut states = na::DMatrix::<f64>::zeros(x0.nrows(), 1);
        states.column_mut(0).copy_from(&x0);

        // // We store the computed inputs
        let inputs = na::DMatrix::<f64>::zeros(0, 0);

        // // # we store the output vectors of the controlled state trajectory
        let mut outputs = na::DMatrix::<f64>::zeros(1, model.get_mat_c().nrows());
        outputs.row_mut(0).copy_from(&(model.get_mat_c() * x0));

        Ok(Controller {
            model: model,
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
