// use mpc_rs::mpc::controller::Controller;
// use mpc_rs::mpc::simulator::system_simulate;

// use plotters::prelude::*;

extern crate nalgebra as na;

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    let mat_ac = na::dmatrix![
        0.0, 1.0, 0.0, 0.0;
        -(k1 + k2) / m1, -(d1 + d2) / m1, k2 / m1, d2 / m1;
        0.0, 0.0, 0.0, 1.0;
        k2 / m2, d2 / m2, -k2 / m2, -d2 / m2
    ];
    let mat_bc = na::dmatrix![0.0; 0.0; 0.0; 1.0 / m2];
    let mat_cc = na::dmatrix![1.0, 0.0, 0.0, 0.0];

    let r = 1usize;
    let m = 1usize; // number of inputs and outputs
    let _n = 4usize; // state dimension

    // Discretization constant
    let sampling = 0.05f64;

    // Model discretization
    let mat_i = na::DMatrix::<f64>::identity(mat_ac.nrows(), mat_ac.nrows());
    let mat_a = mat_ac.scale(sampling).try_inverse().unwrap();
    let mat_b = mat_a * mat_bc.scale(sampling);
    let mat_c = mat_cc;

    // check the eigenvalues
    // let _eigen_a = mat_ac.eig()?;
    // let _eigen_aid = mat_a.eig()?;

    // let time_sample_test = 200;

    // // Compute the system's step response
    // let input_test = 10.0 * Array2::ones((1, time_sample_test));
    // let x0_test = Array::zeros(4);

    // // // # simulate the discrete-time system
    // let (y_test, _x_test) = system_simulate(&mat_a, &mat_b, &mat_c, &input_test, &x0_test);

    // // Draw the response
    // {
    //     let root = BitMapBackend::new("step_response.png", (800, 600)).into_drawing_area();
    //     root.fill(&WHITE)?;

    //     let max_y = y_test.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    //     let min_y = y_test.iter().cloned().fold(f64::INFINITY, f64::min);
    //     let mut chart = ChartBuilder::on(&root)
    //         .caption("System Output Y", ("sans-serif", 20))
    //         .margin(10)
    //         .x_label_area_size(30)
    //         .y_label_area_size(40)
    //         .build_cartesian_2d(0..y_test.ncols() as i32, min_y..max_y)?;

    //     chart.configure_mesh().draw()?;

    //     // Plot input
    //     let series_input: Vec<(i32, f64)> = input_test
    //         .row(0)
    //         .iter()
    //         .enumerate()
    //         .map(|(i, &val)| (i as i32, val as f64))
    //         .collect();

    //     chart
    //         .draw_series(LineSeries::new(series_input, &Palette99::pick(0)))?
    //         .label(format!("Output {}", 0))
    //         .legend(move |(x, y)| PathElement::new([(x, y), (x + 20, y)], &Palette99::pick(0)));

    //     // Plot system response
    //     let series_y: Vec<(i32, f64)> = y_test
    //         .row(0)
    //         .iter()
    //         .enumerate()
    //         .map(|(i, &val)| (i as i32, val as f64))
    //         .collect();

    //     chart.draw_series(LineSeries::new(series_y, &Palette99::pick(1)))?;

    //     chart
    //         .configure_series_labels()
    //         .background_style(&WHITE)
    //         .border_style(&BLACK)
    //         .draw()?;
    // }

    // W1 matrix
    // let mut mat_w1: Array2<f64> = Array2::zeros((v * m, v * m));

    // mat_w1.slice_mut(s![0..m, 0..m]).assign(&Array2::eye(m));

    // for i in 1..v {
    //     mat_w1
    //         .slice_mut(s![i * m..(i + 1) * m, i * m..(i + 1) * m])
    //         .assign(&Array2::eye(m));
    //     mat_w1
    //         .slice_mut(s![i * m..(i + 1) * m, (i - 1) * m..(i) * m])
    //         .assign(&(-1.0 * Array2::eye(m)));
    // }

    // // W2 matrix
    // let mat_q0 = array![0.0000000011f64];
    // let math_q_other = array![0.0001f64];

    // let mut mat_w2: Array2<f64> = Array2::zeros((v * m, v * m));

    // mat_w2.slice_mut(s![0..m, 0..m]).assign(&mat_q0);

    // for i in 1..v {
    //     mat_w2
    //         .slice_mut(s![i * m..(i + 1) * m, i * m..(i + 1) * m])
    //         .assign(&math_q_other);
    // }

    // // W3 matrix
    // let mat_w3 = mat_w1.t().dot(&(mat_w2.dot(&mat_w1)));

    // // W4 matrix
    // let mut mat_w4: Array2<f64> = Array2::zeros((f * r, f * r));

    // let pred_weight = array![10f64];

    // for i in 0..f {
    //     mat_w4
    //         .slice_mut(s![i * r..(i + 1) * r, i * r..(i + 1) * r])
    //         .assign(&pred_weight);
    // }

    // let time_steps = 300;

    // // Define a step trajectory
    // // let desired_traj: Array2<f64> = 0.3 * Array2::ones((time_steps, 1));

    // // Define a pulse trajectory
    // let mut desired_traj: Array2<f64> = Array2::zeros((time_steps, 1));
    // desired_traj
    //     .slice_mut(s![0..100, ..])
    //     .assign(&Array2::ones((100, 1)));
    // desired_traj
    //     .slice_mut(s![200.., ..])
    //     .assign(&Array2::ones((100, 1)));

    // // Set the initial state
    // let x0 = x0_test;

    // // Create the controller
    // let mut mpc = Controller::new(
    //     mat_a,
    //     mat_b,
    //     mat_c,
    //     f,
    //     v,
    //     &mat_w3,
    //     &mat_w4,
    //     x0,
    //     &desired_traj,
    // )?;

    // for _ in 0..time_steps - f {
    //     mpc.compute_control_inputs();
    // }

    // {
    //     let root = BitMapBackend::new("control.png", (800, 600)).into_drawing_area();
    //     root.fill(&WHITE)?;

    //     let max_y = mpc
    //         .outputs
    //         .iter()
    //         .cloned()
    //         .fold(f64::NEG_INFINITY, f64::max);
    //     let min_y = mpc.outputs.iter().cloned().fold(f64::INFINITY, f64::min);
    //     let mut chart = ChartBuilder::on(&root)
    //         .caption("System Output Y", ("sans-serif", 20))
    //         .margin(10)
    //         .x_label_area_size(30)
    //         .y_label_area_size(40)
    //         .build_cartesian_2d(0..y_test.ncols() as i32, min_y..max_y)?;

    //     chart.configure_mesh().draw()?;

    //     // Plot input
    //     let inputs_series: Vec<(i32, f64)> = mpc
    //         .inputs
    //         .row(0)
    //         .iter()
    //         .enumerate()
    //         .map(|(i, &val)| (i as i32, val as f64))
    //         .collect();

    //     chart
    //         .draw_series(LineSeries::new(inputs_series, &Palette99::pick(0)))?
    //         .label(format!("Input"))
    //         .legend(move |(x, y)| PathElement::new([(x, y), (x + 20, y)], &Palette99::pick(0)));

    //     // Plot system response
    //     let outputs_serie: Vec<(i32, f64)> = mpc
    //         .outputs
    //         .iter()
    //         .enumerate()
    //         .map(|(i, &val)| (i as i32, val as f64))
    //         .collect();

    //     chart
    //         .draw_series(LineSeries::new(outputs_serie, &Palette99::pick(1)))?
    //         .label(format!("Output"))
    //         .legend(move |(x, y)| PathElement::new([(x, y), (x + 40, y)], &Palette99::pick(1)));

    //     let traj_series: Vec<(i32, f64)> = mpc
    //         .desired_ctrl_traj_total
    //         .iter()
    //         .enumerate()
    //         .map(|(i, &val)| (i as i32, val as f64))
    //         .collect();

    //     chart
    //         .draw_series(LineSeries::new(traj_series, &Palette99::pick(2)))?
    //         .label(format!("Desired trajectory"))
    //         .legend(move |(x, y)| PathElement::new([(x, y), (x + 60, y)], &Palette99::pick(2)));

    //     chart
    //         .configure_series_labels()
    //         .background_style(&WHITE)
    //         .border_style(&BLACK)
    //         .draw()?;
    // }

    Ok(())
}
