use bevy::{input::mouse::AccumulatedMouseScroll, prelude::*};

use crate::fluid_simulation::{Particle, SimulationParameters};

const BOUNDING_BOX_SIZE: f32 = 25.0;
const BOUNDING_BOX_MAX_SIZE: f32 = 200.0;
const BOUNDING_BOX_MIN_SIZE: f32 = 0.0;
const BOUNDING_BOX_DAMPENING_FACTOR: f32 = 0.8;

pub struct BoundingBox;

impl Plugin for BoundingBox {
    fn build(&self, app: &mut App) {
        app
        .add_systems(PreStartup, spawn_bounding_box)
        .add_systems(Update, (render_bounding_box, control_bounding_box, check_bounding_box_collision));
    }
}

#[derive(Component)]
struct BoundingBoxStruct {
    size: f32,
    min_size: f32,
    max_size: f32
}

// Spawns the bounding box to limit where the particles can go
fn spawn_bounding_box(
    mut commands: Commands
) {
    commands.spawn(BoundingBoxStruct {size: BOUNDING_BOX_SIZE, min_size: BOUNDING_BOX_MIN_SIZE, max_size: BOUNDING_BOX_MAX_SIZE});
}

// Render Bounding Box using gizmos
fn render_bounding_box(
    mut gizmos: Gizmos,
    query: Query<&BoundingBoxStruct>
) {

    let bounding_box = query.single().unwrap();

    gizmos.cuboid(
        Transform::from_xyz(0.0, 0.0, 0.0).with_scale(Vec3::splat(bounding_box.size)),
        Color::WHITE,
    );
}

// Control the size of the bounding box with keys
fn control_bounding_box(
    keys: Res<ButtonInput<KeyCode>>,
    mut query: Query<&mut BoundingBoxStruct>,
    mouse_wheel_input: Res<AccumulatedMouseScroll>,
    keyboard_input: Res<ButtonInput<KeyCode>>,
) {

    let mut bounding_box = query.single_mut().unwrap();

    if keyboard_input.pressed(KeyCode::KeyB) {
        bounding_box.size += mouse_wheel_input.delta.y * 2.0;

    }

    if keyboard_input.pressed(KeyCode::KeyB) && keyboard_input.pressed(KeyCode::ShiftLeft) {
        bounding_box.size += mouse_wheel_input.delta.y * 4.0;

    }
}

// Check the collisions between the particles and the bounding box limits
fn check_bounding_box_collision(
    mut particles: Query<(&mut Transform, &mut Particle)>,
    bounding_box: Query<&BoundingBoxStruct>,
    simulation_parameters: Res<SimulationParameters>,
    time: Res<Time>
) {
    let bounding_box = bounding_box.single().unwrap();

    let particle_size = 0.1;

    let half_bounds_size = (bounding_box.size / 2.0) as f32 - Vec3::new(1.0, 1.0, 1.0) * particle_size;

    for (mut particle_transform, mut particle) in &mut particles {
        if particle_transform.translation.x.abs() > half_bounds_size.x {
            particle_transform.translation.x = half_bounds_size.x * particle_transform.translation.x.signum();
            particle.velocity *= -1.0 * BOUNDING_BOX_DAMPENING_FACTOR;
            particle_transform.translation += particle.velocity * 1.0 / simulation_parameters.time_step;
        }

        if particle_transform.translation.y.abs() > half_bounds_size.y {
            particle_transform.translation.y = half_bounds_size.y * particle_transform.translation.y.signum();
            particle.velocity *= -1.0 * BOUNDING_BOX_DAMPENING_FACTOR;
            particle_transform.translation += particle.velocity * 1.0 / simulation_parameters.time_step;
        }

        if particle_transform.translation.z.abs() > half_bounds_size.z {
            particle_transform.translation.z = half_bounds_size.z * particle_transform.translation.z.signum();
            particle.velocity *= -1.0 * BOUNDING_BOX_DAMPENING_FACTOR;
            particle_transform.translation += particle.velocity * 1.0 / simulation_parameters.time_step;
        }
    }
}