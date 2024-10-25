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

fn main() {
    let mut system = simulator::TwoWheeledSystem{x: 0.0, y: 0.0, yaw: 0.0, v: 0.0};

    let command = simulator::Command{u1: 1.0, u2: 0.0};
    let delta_t = 0.1;

    println!("Initial state : {:?}", system);

    for _ in 1..100
    {
        simulator::apply_differential_model(&mut system, &command, delta_t);
        println!("{:?}", system);
    }
}

