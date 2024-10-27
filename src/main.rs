use plotters::prelude::*;

mod simulator
{
    #[derive(Debug, Copy, Clone)]
    pub struct TwoWheeledSystem
    {
        pub x: f32,
        pub y: f32,
        pub yaw: f32,
        pub v: f32
    }

    pub struct Command
    {
        pub u1: f32,
        pub u2: f32
    }

    pub fn apply_differential_model(system: &mut TwoWheeledSystem, cmd: &Command, dt: f32)
    {
        let dx = system.yaw.cos() * system.v;
        let dy = system.yaw.sin() * system.v;
        let dv = cmd.u1;
        let d_yaw = system.v / 0.25 * cmd.u2.sin();

        system.x += dt * dx;
        system.y += dt * dy;
        system.yaw += dt * d_yaw;
        system.v += dt * dv;
    }
}

use ndarray::prelude::*;
use ndarray::Array;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define parameters
    let m1 = 2.0;
    let m2 = 2.0;
    let k1 = 100.0;
    let k2 = 200.0; 
    let d1 = 1.0;
    let d2 = 5.0; 
    // Define the continuous-time system matrices
    let ac = array![[0.0, 1.0, 0.0, 0.0],
                    [-(k1+k2)/m1 ,  -(d1+d2)/m1 , k2/m1 , d2/m1 ],
                    [0.0 , 0.0 ,  0.0 , 1.0], 
                    [k2/m2,  d2/m2, -k2/m2, -d2/m2]];
    let bc = array![[0.0],[0.0],[0.0],[1.0/m2]];
    let cc = array![[1.0, 0.0, 0.0, 0.0]];

    let r = 1;
    let m = 1; // number of inputs and outputs
    let n = 4; // state dimension

    // Discretization constant
    let sampling = 0.05;

    // Model discretization
    let I : Array::<f64, _> = Array::eye(ac.shape()[0]);
    A=np.linalg.inv(I-sampling*Ac)
    B=A*sampling*Bc
    C=Cc

    // # check the eigenvalues
    // eigen_A=np.linalg.eig(Ac)[0]
    // eigen_Aid=np.linalg.eig(A)[0]

    // timeSampleTest=200

    // # compute the system's step response
    // inputTest=10*np.ones((1,timeSampleTest))
    // x0test=np.zeros(shape=(4,1))


    // # simulate the discrete-time system 
    // Ytest, Xtest=systemSimulate(A,B,C,inputTest,x0test)

    let mut system = simulator::TwoWheeledSystem{x: 0.0, y: 0.0, yaw: 0.0, v: 0.0};

    let command = simulator::Command{u1: 1.0, u2: 0.0};
    let delta_t = 0.1;

    println!("Initial state : {:?}", system);

    let mut trajectory = Vec::new();
    for _ in 1..100
    {
        simulator::apply_differential_model(&mut system, &command, delta_t);
        trajectory.push(system.clone());
        println!("{:?}", system);
    }

    let root = BitMapBackend::new("./test.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-50f32..50f32, -0.1f32..1f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            trajectory.iter().map(|state| (state.x, state.y)),
            &RED
        ))?;

    root.present()?;

    Ok(())
}

