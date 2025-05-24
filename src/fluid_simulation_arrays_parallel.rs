use rand::{seq::index, Rng};
use std::{cell, f32::consts::PI};

use bevy::{
    color::palettes::css::{BLUE, RED},
    math::VectorSpace,
    prelude::*, tasks::{ComputeTaskPool, ParallelSliceMut},
};

const GRAVITY: f32 = 9.81;

const TARGET_DENSITY: f32 = 1.0;

const PRESSURE_MULTIPLIER: f32 = 1.0;

const PARTICLE_NUMBER: u32 = 4000;
const PARTICLE_SPACING: f32 = 0.1;
const PARTICLE_SIZE: f32 = 0.1;
const PARTICLE_INITIAL_VELOCITY: Vec3 = Vec3::new(0.0, 0.0, 0.0);
const PARTICLE_INITIAL_DENSITY: f32 = TARGET_DENSITY;
const PARTICLE_MASS: f32 = 1.0;
const PARTICLE_SMOOTHING_RADIUS: f32 = 1.35;

pub struct FluidSimulation;

impl Plugin for FluidSimulation {
    fn build(&self, app: &mut bevy::app::App) {
        app.insert_resource(SimulationParameters {
            gravity: GRAVITY,
            target_density: TARGET_DENSITY,
            pressure_multiplier: PRESSURE_MULTIPLIER,
            particle_smoothing_radius: PARTICLE_SMOOTHING_RADIUS,
        })
        .add_systems(Startup, spawn_particles)
        .add_systems(Update, simulate_fluid)
        .add_systems(Update, visual_debug);
    }
}

#[derive(Component)]
pub struct Particle {
    pub id: u32,
    pub velocity: Vec3,
    pub mass: f32,
    pub density: f32,
    pub pressure_force: Vec3,
    pub position: Vec3,
    pub cell_key: u32
}

#[derive(Component)]
pub struct SmoothingRadiusDebug;

#[derive(Resource)]
pub struct SimulationParameters {
    pub gravity: f32,
    pub target_density: f32,
    pub pressure_multiplier: f32,
    pub particle_smoothing_radius: f32,
}

// Spawns all the fluid particles
fn spawn_particles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    let particles_per_row: f32 = f32::sqrt(PARTICLE_NUMBER as f32);
    let particles_per_col: f32 = (PARTICLE_NUMBER - 1) as f32 / particles_per_row + 1.0;
    let spacing: f32 = PARTICLE_SIZE * 2.0 + PARTICLE_SPACING;

    let mut x: f32;
    let mut y: f32;

    for i in 0..PARTICLE_NUMBER {
        x = (i as f32 % particles_per_row - particles_per_row / 2.0).round() * spacing;
        y = (i as f32 / particles_per_row - particles_per_col / 2.0).round() * spacing;

        //println!("x: {x} y: {y}");

        commands.spawn((
            Particle {
                id: i,
                velocity: PARTICLE_INITIAL_VELOCITY,
                mass: PARTICLE_MASS,
                density: PARTICLE_INITIAL_DENSITY,
                pressure_force: Vec3::ZERO,
                position: Vec3::new(rand::rng().random_range(-12.0..=12.0), rand::rng().random_range(-12.0..=12.0), 0.0),
                cell_key: 0
            },
            Mesh2d(meshes.add(Circle::new(PARTICLE_SIZE))),
            MeshMaterial2d(materials.add(ColorMaterial::from_color(RED))),
            Transform::from_xyz(rand::rng().random_range(-12.0..=12.0), rand::rng().random_range(-12.0..=12.0), 0.0),
        ))/* .with_child((
            Transform::from_scale(Vec3::new(0.025, 0.025, 1.0)),
            Text2d::new(""),
        )).with_child((
            Mesh2d(meshes.add(Circle::new(PARTICLE_SMOOTHING_RADIUS))),
            MeshMaterial2d(materials.add(ColorMaterial::from_color(Srgba::new(0.0, 0.0, 1.0, 0.05)))),
            SmoothingRadiusDebug
        )  ) */;
    }
}

fn visual_debug(mut query: Query<(&Particle, &Children)>, mut text_query: Query<&mut Text2d>) {
    for (particle, children) in query.iter_mut() {
        for child in children.iter() {
            if let Ok(mut text) = text_query.get_mut(*child) {

                /* if particle.id % 100 == 0 {
                    text.0 = format!("id: {}\ndensity {:.2}\npressure force {:.2}\nvelocity {:.2}\ncell key {}", particle.id, particle.density, particle.pressure_force, particle.velocity, particle.cell_key);
                } */

                // text.0 = format!("{:.2}", particle.density);
                // text.0 = format!("{:.2}", particle.pressure_force);
                //text.0 = format!("{:.2}", particle.id);
                // text.0 = format!("{:.2}", particle.velocity);
                //text.0 = format!("{:.2}", particle.position);
                // text.0 = format!("{:.2}", particle.cell_key);
            }
        }
    }
}

fn smoothing_kernel(smoothing_radius: f32, distance: f32) -> f32 {
    if distance >= smoothing_radius  {
        return 0.0;
    }

    let volume = PI * f32::powf(smoothing_radius, 4.0) / 6.0;
    return (smoothing_radius - distance) * (smoothing_radius - distance) / volume;
}

fn smoothing_kernel_derivative(smoothing_radius: f32, distance: f32) -> f32 {
    if distance >= smoothing_radius {
        return 0.0;
    }

    let scale: f32 = 12.0 / (PI * f32::powf(smoothing_radius, 4.0));
    return (distance - smoothing_radius) * scale;
}

// Converts a density value to a pressure value
fn convert_density_to_pressure(
    density: f32,
    simulation_parameters: &Res<SimulationParameters>,
) -> f32 {
    let density_error: f32 = density - simulation_parameters.target_density;
    // DEBUG
    //println!("density: {} target_density: {}", density, simulation_parameters.target_density);
    let pressure: f32 = density_error * simulation_parameters.pressure_multiplier;
    return pressure;
}

// Simulate the fluid
fn simulate_fluid(
    mut particles_transform: Query<&mut Transform, With<Particle>>,
    mut particles: Query<&mut Particle>,
    time: Res<Time>,
    simulation_parameters: Res<SimulationParameters>,
) {
    // Create all the buffers that will serve as a cache for the particle values
    let mut cached_positions: [Vec3; PARTICLE_NUMBER as usize] = [Vec3::ZERO; PARTICLE_NUMBER as usize];

    let mut cached_predicted_positions: [Vec3; PARTICLE_NUMBER as usize] = [Vec3::ZERO; PARTICLE_NUMBER as usize];

    let mut cached_velocities: [Vec3; PARTICLE_NUMBER as usize] = [Vec3::ZERO; PARTICLE_NUMBER as usize];

    let mut cached_densities: [f32; PARTICLE_NUMBER as usize] = [TARGET_DENSITY; PARTICLE_NUMBER as usize];

    let mut cached_pressure_forces: [Vec3; PARTICLE_NUMBER as usize] = [Vec3::ZERO; PARTICLE_NUMBER as usize];

    let mut cached_cell_coords: [(i32, i32); PARTICLE_NUMBER as usize] = [(0, 0); PARTICLE_NUMBER as usize];

    let mut cached_cell_keys: [u32; PARTICLE_NUMBER as usize] = [u32::MAX; PARTICLE_NUMBER as usize];

    let mut spatial_lookup: [(usize, u32); PARTICLE_NUMBER as usize] = [(0, 0); PARTICLE_NUMBER as usize];

    let mut start_indices: [u32; PARTICLE_NUMBER as usize] = [u32::MAX; PARTICLE_NUMBER as usize];

    let cell_offsets: [(i32, i32); 9] = [(-1, 1), (0, 1), (1, 1), (-1, 0), (0, 0), (1, 0), (-1, -1), (0, -1), (1, -1)];

    // Loop over all particles to save their positions on the cache buffer
    for (i, particle_transform) in particles_transform.iter().enumerate() {

        // Save all the particle positions
        cached_positions[i] = particle_transform.translation;

        // You can use predicted positions instead to improve simulation results
        cached_predicted_positions[i] = cached_positions[i] + cached_velocities[i] * time.delta_secs();

        // Update spatial lookup buffer
        // Calculate the cell key for each particle based on its current position

        // Get the cell coordinates
        let cell_x = (cached_predicted_positions[i].x / simulation_parameters.particle_smoothing_radius) as i32;
        let cell_y = (cached_predicted_positions[i].y / simulation_parameters.particle_smoothing_radius) as i32;

        // Save the cell coordinates to the cache buffer
        cached_cell_coords[i] = (cell_x, cell_y);

        // Hash the cell coordinates to a single number (this makes easier to work with)
        let cell_hash_x = cell_x * 15823;
        let cell_hash_y = cell_y * 9737333;
        let cell_hash = (cell_hash_x + cell_hash_y) as u32;

        // Get the cell key
        let cell_key = cell_hash % PARTICLE_NUMBER;

        // Save the cell key for each particle into the cache buffer
        cached_cell_keys[i] = cell_key;

        // Save the particle index and its current cell key to the cache buffer
        spatial_lookup[i] = (i, cell_key);

    }

    // Sort the spatial lookup array ordered by the cell key
    spatial_lookup.sort_by_key(|k| k.1);

    // Calculate the start indice of each cell (the starting position on the spatial lookup array from where the particles are stored for each cell)
    for i in 0..PARTICLE_NUMBER as usize {

        let key = spatial_lookup[i].1;

        let key_prev;

        if i == 0 {
            key_prev = u32::MAX;
        } else {
            key_prev = spatial_lookup[i - 1].1;
        }

        if key != key_prev {
            start_indices[key as usize] = i as u32;
        }

    }

    cached_densities.par_splat_map_mut(ComputeTaskPool::get(), None, |i, densities| {
        for (j, density) in densities.iter_mut().enumerate() {

            let mut new_density = 1.0;

        // Get the center of the 3x3 grid cell keys bounding box (inside of which we need to check for the influence of other particles)
        let center_x = cached_cell_coords[j].0;
        let center_y = cached_cell_coords[j].1;

        let sqr_radius = simulation_parameters.particle_smoothing_radius * simulation_parameters.particle_smoothing_radius;

        for (cell_offset_x, cell_offset_y) in cell_offsets {

            // Get the cell hash of the current offset grid cell
            let cell_hash_x = (center_x + cell_offset_x) * 15823;
            let cell_hash_y = (center_y + cell_offset_y) * 9737333;
            let cell_hash = (cell_hash_x + cell_hash_y) as u32;

            // Get the cell key of the current offset grid cell
            let cell_key = cell_hash % PARTICLE_NUMBER;

            // Get the starting position of the spatial lookup array from where the particles inside this offset grid cell are stored
            let cell_start_index = start_indices[cell_key as usize] as usize;

            // Check all the particles inside the current offset grid cell
            for k in cell_start_index..PARTICLE_NUMBER as usize {

                // Break loop if the cell key is different (if the cell key is different it means that the particle is outside this offset grid cell)
                if spatial_lookup[k].1 != cell_key {
                    break;
                }

                // Get the particle index
                let particle_index = spatial_lookup[k].0;

                let sqr_dst = Vec3::length_squared(cached_predicted_positions[particle_index] - cached_predicted_positions[j]);

                // Test if the point is inside the radius
                if sqr_dst <= sqr_radius {

                    let distance = Vec3::length( cached_predicted_positions[particle_index] - cached_predicted_positions[j]);

                    let influence = smoothing_kernel(simulation_parameters.particle_smoothing_radius, distance);

                    new_density = new_density + PARTICLE_MASS * influence;
                    
                }
            }

        }

        *density = new_density;

        }
    });

    cached_pressure_forces.par_splat_map_mut(ComputeTaskPool::get(), None, |i, pressure_forces| {
        for (j, pressure_force) in pressure_forces.iter_mut().enumerate() {

            let mut new_pressure_force = Vec3::ZERO;

        // Get the center of the 3x3 grid cell keys bounding box (inside of which we need to check for the influence of other particles)
        let center_x = cached_cell_coords[j].0;
        let center_y = cached_cell_coords[j].1;

        let sqr_radius = simulation_parameters.particle_smoothing_radius * simulation_parameters.particle_smoothing_radius;

        for (cell_offset_x, cell_offset_y) in cell_offsets {

            // Get the cell hash of the current offset grid cell
            let cell_hash_x = (center_x + cell_offset_x) * 15823;
            let cell_hash_y = (center_y + cell_offset_y) * 9737333;
            let cell_hash = (cell_hash_x + cell_hash_y) as u32;

            // Get the cell key of the current offset grid cell
            let cell_key = cell_hash % PARTICLE_NUMBER;

            // Get the starting position of the spatial lookup array from where the particles inside this offset grid cell are stored
            let cell_start_index = start_indices[cell_key as usize] as usize;

            // Check all the particles inside the current offset grid cell
            for k in cell_start_index..PARTICLE_NUMBER as usize {

                // Break loop if the cell key is different (if the cell key is different it means that the particle is outside this offset grid cell)
                if spatial_lookup[k].1 != cell_key {
                    break;
                }

                // Get the particle index
                let particle_index = spatial_lookup[k].0;

                // Skip calculating pressure on itself
                if  i == particle_index {continue;}

                let sqr_dst = Vec3::length_squared(cached_predicted_positions[particle_index] - cached_predicted_positions[j]);

                // Test if the point is inside the radius
                if sqr_dst <= sqr_radius {

                    let distance = Vec3::length( cached_predicted_positions[particle_index] - cached_predicted_positions[j]);

                    let direction;

                    // If 2 particles are on the exact same spot choose a random direction to apply the pressure force
                    if distance == 0.0 {

                        direction = Vec3::new(
                            rand::rng().random_range(-1.0..=1.0),
                            rand::rng().random_range(-1.0..=1.0),
                            0.0,
                        );
                        
                    } else {

                        direction = (cached_predicted_positions[particle_index] - cached_predicted_positions[j]) / distance;

                    }

                    let slope = smoothing_kernel_derivative(simulation_parameters.particle_smoothing_radius, distance);

                    let new_pressure_force_a = convert_density_to_pressure(cached_densities[j], &simulation_parameters);
                    let new_pressure_force_b = convert_density_to_pressure(cached_densities[particle_index], &simulation_parameters);

                    let shared_new_pressure_force = (new_pressure_force_a + new_pressure_force_b) / 2.0;

                    new_pressure_force += shared_new_pressure_force * direction * slope * PARTICLE_MASS / cached_densities[particle_index];

                }
            }

        }

        *pressure_force = new_pressure_force;

        }});

    // Apply the calculated pressure force to all the particles
    for (i, mut particle_transform) in particles_transform.iter_mut().enumerate() {

        // Store each particle velocity in the cache buffer
        cached_velocities[i] += (cached_pressure_forces[i] / cached_densities[i]) * time.delta_secs();

        // Aply gravity
        // cached_velocities[i] += Vec3::NEG_Y * simulation_parameters.gravity * time.delta_secs();
        
        // Apply movement to the bevy transform
        particle_transform.translation += cached_velocities[i] * time.delta_secs();

    }

    // Set all the debug info for each particle
    for (i, mut particle) in particles.iter_mut().enumerate() {

        particle.position = cached_positions[i];
        particle.density = cached_densities[i];
        particle.velocity = cached_velocities[i];
        particle.pressure_force = cached_pressure_forces[i];
        particle.cell_key = cached_cell_keys[i];


    }


}
