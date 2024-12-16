use mpc_rs::mpc::controller::Controller;
use mpc_rs::mpc::example::generate_pulse_trajectory; //, Model};
use mpc_rs::mpc::dc_motor_model::Model;
use mpc_rs::mpc::simulator::compute_system_response;

use plotters::prelude::*;

use std::borrow::Borrow;
use std::rc::Rc;

extern crate nalgebra as na;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let J = 0.01f64; // [kg.m^2] Moment of inertia of the robot
    let b = 0.1f64;  // [N.m.s] Motor discous friction constant
    let Ke = 0.01f64; // [V/rad/sec] Electromotive force constant
    let Kt = 0.01f64; // [N.m/Amp] Motor torque constant
    let R = 1.0f64; // [Ohm] Electric resistance       
    let L = 0.5f64; // [H] Electric inductance 

    let sampling_dt = 0.05;

    let model = Rc::new(Model::new(b, J, Ke, L, R, sampling_dt));

    let sampling_time = 10.0f64;
    let n_samples = (sampling_time / sampling_dt).floor() as usize;
    let input_test = na::DMatrix::from_element(1, n_samples, 10.0f64);
    let x0_test = na::DVector::<f64>::zeros(2);

    let system_response = compute_system_response(
        <Rc<Model> as Borrow<Model>>::borrow(&model),
        &input_test,
        &x0_test,
    );

    // Define parameters
    let f = 20usize;
    let v = 20usize;
    let time_steps = 300;

    // Define a pulse trajectory
    let trajectory = generate_pulse_trajectory(time_steps);

    // // Set the initial state
    // let x0 = x0_test;

    // // Create the controller
    // let pred_weight = 10f64;
    // let q0 = 0.0000000011f64;
    // let q_other = 0.0001f64;

    // let mut mpc = Controller::new(model, f, v, q0, q_other, pred_weight, x0, &trajectory)?;

    // for _ in 0..time_steps - f {
    //     mpc.compute_control_inputs();
    // }

    // Draw the step response
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

    // Draw the control
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
    //         .build_cartesian_2d(0..system_response.len() as i32, min_y..max_y)?;

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
