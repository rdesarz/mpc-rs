extern crate nalgebra as na;

pub trait LinearDiscreteModel {
    fn get_mat_a(&self) -> &na::DMatrix<f64>;
    fn get_mat_b(&self) -> &na::DMatrix<f64>;
    fn get_mat_c(&self) -> &na::DMatrix<f64>;
    fn get_sampling_dt(&self) -> f64;
}
