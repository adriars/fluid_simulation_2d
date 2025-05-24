use bevy::{
    input::mouse::{AccumulatedMouseScroll, MouseScrollUnit, MouseWheel},
    prelude::*,
};

use crate::fluid_simulation::{Particle, SimulationParameters, SmoothingRadiusDebug, MAX_TIME_STEP, MIN_TIME_STEP};

pub struct UIDebug;

#[derive(Component)]
struct GravityUiNode;
#[derive(Component)]
struct TargetDensityUiNode;
#[derive(Component)]
struct PressureMultiplierUiNode;
#[derive(Component)]
struct ParticleSmoothingRadiusUiNode;
#[derive(Component)]
struct ParticleTimeStepUiNode;

impl Plugin for UIDebug {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, spawn_ui_debug_nodes)
            .add_systems(Update, update_ui_debug_nodes);
    }
}

fn spawn_ui_debug_nodes(mut commands: Commands,
    simulation_parameters: Res<SimulationParameters>
) {
    commands
        .spawn(
            (Node {
                // fill the entire window
                width: Val::Percent(100.),
                height: Val::Percent(100.),
                flex_direction: FlexDirection::Column,
                align_items: AlignItems::Start,
                ..Default::default()
            }),
        )
        .with_child((
            Text::new(format!("Gravity: {} (G)", simulation_parameters.gravity)),
            GravityUiNode,
        ))
        .with_child((
            Text::new(format!("Target density: {} (D)", simulation_parameters.target_density)),
            TargetDensityUiNode,
        ))
        .with_child((
            Text::new(format!("Pressure multiplier: {} (P)", simulation_parameters.pressure_multiplier)),
            PressureMultiplierUiNode,
        ))
        .with_child((
            Text::new(format!("Particle smoothing radius: {} (S)", simulation_parameters.particle_smoothing_radius)),
            ParticleSmoothingRadiusUiNode,
        )).with_child((
            Text::new(format!("Time step {} (T)", simulation_parameters.time_step)),
            ParticleTimeStepUiNode,
        ))
        .with_child(Text::new("Bounding box (B)"))
        .with_child(Text::new("Zoom (Z)"))
        .with_child(
            (Node {
                // fill the entire window
                width: Val::Percent(100.),
                height: Val::Percent(100.),
                flex_direction: FlexDirection::Column,
                align_items: AlignItems::End,
                ..Default::default()
            })
        ).with_child(Text::new("Hold the legend key + mouse scroll wheel to modify the parameters (hold left shift to bigger steps)"));

}

fn update_gravity_ui_node(
    gravity_ui_node: &mut Query<&mut Text, With<GravityUiNode>>,
    simulation_parameters: &ResMut<SimulationParameters>
) {
    gravity_ui_node.single_mut().unwrap().0 = format!(
        "Gravity: {} (G)",
        simulation_parameters.gravity
    );
}

fn update_target_density_ui_node(
    target_density_ui_node: &mut Query<&mut Text, With<TargetDensityUiNode>>,
    simulation_parameters: &ResMut<SimulationParameters>
) {
    target_density_ui_node.single_mut().unwrap().0 = format!(
        "Target density: {} (D)",
        simulation_parameters.target_density
    );
}

fn update_pressure_multiplier_ui_node(
    pressure_multiplier_ui_node: &mut Query<&mut Text, With<PressureMultiplierUiNode>>,
    simulation_parameters: &ResMut<SimulationParameters>
) {
    pressure_multiplier_ui_node.single_mut().unwrap().0 = format!(
        "Pressure multiplier: {} (P)",
        simulation_parameters.pressure_multiplier
    );
}

fn update_particle_smoothing_radius_ui_node(
    particle_smoothing_radius_ui_node: &mut Query<&mut Text, With<ParticleSmoothingRadiusUiNode>>,
    simulation_parameters: &ResMut<SimulationParameters>
) {
    particle_smoothing_radius_ui_node.single_mut().unwrap().0 = format!(
        "Particle smoothing radius: {} (S)",
        simulation_parameters.particle_smoothing_radius
    );
}

fn update_time_step_ui_node(
    particle_time_step_ui_node: &mut Query<&mut Text, With<ParticleTimeStepUiNode>>,
    simulation_parameters: &ResMut<SimulationParameters>
) {
    particle_time_step_ui_node.single_mut().unwrap().0 = format!(
        "Time step {} (T)",
        simulation_parameters.time_step
    );
}

fn update_ui_debug_nodes(
    mouse_wheel_input: Res<AccumulatedMouseScroll>,
    keyboard_input: Res<ButtonInput<KeyCode>>,
    mut simulation_parameters: ResMut<SimulationParameters>,
    mut ui_nodes: ParamSet<(
        Query<&mut Text, With<GravityUiNode>>,
        Query<&mut Text, With<TargetDensityUiNode>>,
        Query<&mut Text, With<PressureMultiplierUiNode>>,
        Query<&mut Text, With<ParticleSmoothingRadiusUiNode>>,
        Query<&mut Text, With<ParticleTimeStepUiNode>>
    )>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<&mut Mesh2d, With<SmoothingRadiusDebug>>
) {

    // Gravity UI node
    if keyboard_input.pressed(KeyCode::KeyG) && keyboard_input.pressed(KeyCode::ShiftLeft) {
        simulation_parameters.gravity += mouse_wheel_input.delta.y * 1.0;

        update_gravity_ui_node(&mut ui_nodes.p0(), &simulation_parameters);
    }

    if keyboard_input.pressed(KeyCode::KeyG) {
        simulation_parameters.gravity += mouse_wheel_input.delta.y * 0.1;

        update_gravity_ui_node(&mut ui_nodes.p0(), &simulation_parameters);
    }

    // Target density UI node
    if keyboard_input.pressed(KeyCode::KeyD) && keyboard_input.pressed(KeyCode::ShiftLeft) {
        simulation_parameters.target_density += mouse_wheel_input.delta.y * 1.0;

        update_target_density_ui_node(&mut ui_nodes.p1(), &simulation_parameters);
    }

    if keyboard_input.pressed(KeyCode::KeyD) {
        simulation_parameters.target_density += mouse_wheel_input.delta.y * 0.1;

        update_target_density_ui_node(&mut ui_nodes.p1(), &simulation_parameters);
    }

    // Pressure multiplier UI node
    if keyboard_input.pressed(KeyCode::KeyP) && keyboard_input.pressed(KeyCode::ShiftLeft) {
        simulation_parameters.pressure_multiplier += mouse_wheel_input.delta.y * 1.0;

        update_pressure_multiplier_ui_node(&mut ui_nodes.p2(), &simulation_parameters);
    }

    if keyboard_input.pressed(KeyCode::KeyP) {
        simulation_parameters.pressure_multiplier += mouse_wheel_input.delta.y * 0.1;

        update_pressure_multiplier_ui_node(&mut ui_nodes.p2(), &simulation_parameters);
    }

    // Particle smoothing radius UI node
    if keyboard_input.pressed(KeyCode::KeyS) && keyboard_input.pressed(KeyCode::ShiftLeft) {
        simulation_parameters.particle_smoothing_radius += mouse_wheel_input.delta.y * 1.0;

        update_particle_smoothing_radius_ui_node(&mut ui_nodes.p3(), &simulation_parameters);

       for mut mesh2d in query.iter_mut() {
            mesh2d.0 = meshes.add(Circle::new(simulation_parameters.particle_smoothing_radius));
       }
    }

    if keyboard_input.pressed(KeyCode::KeyS) {
        simulation_parameters.particle_smoothing_radius += mouse_wheel_input.delta.y * 0.1;

        update_particle_smoothing_radius_ui_node(&mut ui_nodes.p3(), &simulation_parameters);

        for mut mesh2d in query.iter_mut() {
            mesh2d.0 = meshes.add(Circle::new(simulation_parameters.particle_smoothing_radius));
       }
    }

    // Simulation speed UI node
    if keyboard_input.pressed(KeyCode::KeyT) && keyboard_input.pressed(KeyCode::ShiftLeft) {

            simulation_parameters.time_step = (simulation_parameters.time_step + mouse_wheel_input.delta.y * 1.0).clamp(MIN_TIME_STEP, MAX_TIME_STEP);

        update_time_step_ui_node(&mut ui_nodes.p4(), &simulation_parameters);

    }

    if keyboard_input.pressed(KeyCode::KeyT) {

            simulation_parameters.time_step = (simulation_parameters.time_step + mouse_wheel_input.delta.y * 0.1).clamp(MIN_TIME_STEP, MAX_TIME_STEP);


        update_time_step_ui_node(&mut ui_nodes.p4(), &simulation_parameters);

    }
}
