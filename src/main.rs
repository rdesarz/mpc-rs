use plotters::prelude::*;

mod simulator
{
    #[derive(Debug)]
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut system = simulator::TwoWheeledSystem{x: 0.0, y: 0.0, yaw: 0.0, v: 0.0};

    let command = simulator::Command{u1: 1.0, u2: 0.0};
    let delta_t = 0.1;

    println!("Initial state : {:?}", system);

    for _ in 1..100
    {
        simulator::apply_differential_model(&mut system, &command, delta_t);
        println!("{:?}", system);
    }

    let root = BitMapBackend::new("./test.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-1f32..1f32, -0.1f32..1f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            &RED,
        ))?
        .label("y = x^2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

