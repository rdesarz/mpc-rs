use mpc_rs::mpc::controller::Controller;
use mpc_rs::mpc::example::Model;
use mpc_rs::mpc::linear_discrete_model::LinearDiscreteModel;
use mpc_rs::mpc::simulator::compute_system_response;

use plotters::prelude::*;

extern crate nalgebra as na;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define parameters
    let f = 20usize;
    let v = 20usize;
    let r = 1usize;
    let m = 1usize; // number of inputs and outputs

    let m1 = 2.0;
    let m2 = 2.0;
    let k1 = 100.0;
    let k2 = 200.0;
    let d1 = 1.0;
    let d2= 5.0; 
    let sampling_dt = 0.05;

    let model = Model::new(m1, m2, k1, k2, d1, d2, sampling_dt);

    let sampling_time = 10.0f64;
    let n_samples = (sampling_time / sampling_dt).floor() as usize;
    let input_test = na::DMatrix::from_element(1, n_samples, 10.0f64);
    let x0_test = na::DVector::<f64>::zeros(4);

    let system_response = compute_system_response(&model, &input_test, &x0_test);

    // Draw the response
    {
        let root = BitMapBackend::new("step_response.png", (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        let max_y = system_response
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_y = system_response
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let mut chart = ChartBuilder::on(&root)
            .caption("System Output Y", ("sans-serif", 20))
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(40)
            .build_cartesian_2d(0..system_response.len() as i32, min_y..max_y)?;

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
        let series_y: Vec<(i32, f64)> = system_response
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
    }

    // W1 matrix
    let mut mat_w1 = na::DMatrix::<f64>::zeros(v * m, v * m);

    mat_w1
        .view_range_mut(0..m, 0..m)
        .copy_from(&na::DMatrix::identity(m, m));

    for i in 1..v {
        mat_w1
            .view_range_mut(i * m..(i + 1) * m, i * m..(i + 1) * m)
            .copy_from(&na::DMatrix::identity(m, m));
        mat_w1
            .view_range_mut(i * m..(i + 1) * m, (i - 1) * m..i * m)
            .copy_from(&(na::DMatrix::identity(m, m).scale(-1.0)));
    }

    // W2 matrix
    let mat_q0 = na::DMatrix::from_element(m, m, 0.0000000011f64);
    let mat_q_other = na::DMatrix::from_element(m, m, 0.0001f64);

    let mut mat_w2 = na::DMatrix::<f64>::zeros(v * m, v * m);

    mat_w2.view_range_mut(0..m, 0..m).copy_from(&mat_q0);

    for i in 1..v {
        mat_w2
            .view_range_mut(i * m..(i + 1) * m, i * m..(i + 1) * m)
            .copy_from(&mat_q_other);
    }

    // W3 matrix
    let mat_w3 = mat_w1.transpose() * (mat_w2 * mat_w1);

    // W4 matrix
    let mut mat_w4 = na::DMatrix::<f64>::zeros(f * r, f * r);

    let pred_weight = 10f64;

    for i in 0..f {
        mat_w4
            .view_range_mut(i * r..(i + 1) * r, i * r..(i + 1) * r)
            .copy_from(&na::DMatrix::from_element(r, r, pred_weight));
    }

    let time_steps = 300;

    // Define a pulse trajectory
    let mut desired_traj = na::DMatrix::<f64>::zeros(time_steps, 1);
    desired_traj
        .view_range_mut(0..100, 0..1)
        .copy_from(&na::DMatrix::from_element(100, 1, 1.0));
    desired_traj
        .view_range_mut(200..300, 0..1)
        .copy_from(&na::DMatrix::from_element(100, 1, 1.0));

    // Set the initial state
    let x0 = x0_test;

    // Create the controller
    let mut mpc = Controller::new(
        model.get_mat_a().clone(),
        model.get_mat_b().clone(),
        model.get_mat_c().clone(),
        f,
        v,
        &mat_w3,
        &mat_w4,
        x0,
        &desired_traj,
    )?;

    for _ in 0..time_steps - f {
        mpc.compute_control_inputs();
    }

    {
        let root = BitMapBackend::new("control.png", (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        let max_y = mpc
            .outputs
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_y = mpc.outputs.iter().cloned().fold(f64::INFINITY, f64::min);
        let mut chart = ChartBuilder::on(&root)
            .caption("System Output Y", ("sans-serif", 20))
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(40)
            .build_cartesian_2d(0..system_response.len() as i32, min_y..max_y)?;

        chart.configure_mesh().draw()?;

        // Plot input
        let inputs_series: Vec<(i32, f64)> = mpc
            .inputs
            .row(0)
            .iter()
            .enumerate()
            .map(|(i, &val)| (i as i32, val as f64))
            .collect();

        chart
            .draw_series(LineSeries::new(inputs_series, &Palette99::pick(0)))?
            .label(format!("Input"))
            .legend(move |(x, y)| PathElement::new([(x, y), (x + 20, y)], &Palette99::pick(0)));

        // Plot system response
        let outputs_serie: Vec<(i32, f64)> = mpc
            .outputs
            .iter()
            .enumerate()
            .map(|(i, &val)| (i as i32, val as f64))
            .collect();

        chart
            .draw_series(LineSeries::new(outputs_serie, &Palette99::pick(1)))?
            .label(format!("Output"))
            .legend(move |(x, y)| PathElement::new([(x, y), (x + 40, y)], &Palette99::pick(1)));

        let traj_series: Vec<(i32, f64)> = mpc
            .desired_ctrl_traj_total
            .iter()
            .enumerate()
            .map(|(i, &val)| (i as i32, val as f64))
            .collect();

        chart
            .draw_series(LineSeries::new(traj_series, &Palette99::pick(2)))?
            .label(format!("Desired trajectory"))
            .legend(move |(x, y)| PathElement::new([(x, y), (x + 60, y)], &Palette99::pick(2)));

        chart
            .configure_series_labels()
            .background_style(&WHITE)
            .border_style(&BLACK)
            .draw()?;
    }

    Ok(())
}
