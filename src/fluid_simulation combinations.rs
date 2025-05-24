use rand::{seq::index, Rng};
use std::{cell, f32::consts::PI};

use bevy::{
    color::palettes::css::{BLUE, RED},
    math::{ops::powf, VectorSpace},
    prelude::*,
};

const GRAVITY: f32 = 9.81;

const TARGET_DENSITY: f32 = 1.0;

const PRESSURE_MULTIPLIER: f32 = 1.0;

const PARTICLE_NUMBER: u32 = 1000;
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
        // .add_systems(Update, simulate_fluid)
        .add_systems(Update, simulate)
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
            Transform::from_scale(Vec3::new(0.1, 0.1, 1.0)),
            Text2d::new(""),
        )) *//* .with_child((
            Mesh2d(meshes.add(Circle::new(PARTICLE_SMOOTHING_RADIUS))),
            MeshMaterial2d(materials.add(ColorMaterial::from_color(Srgba::new(0.0, 0.0, 1.0, 0.05)))),
            SmoothingRadiusDebug
        )  ) */;
    }
}

fn visual_debug(mut query: Query<(&Particle, &Children)>, mut text_query: Query<&mut Text2d>) {
    for (particle, children) in query.iter_mut() {
        for child in children.iter() {
            if let Ok(mut text) = text_query.get_mut(child) {

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

fn simulate (
    mut query: Query<(&mut Particle, &mut Transform)>,
    time: Res<Time>,
    simulation_parameters: Res<SimulationParameters>,
) {

    // Reset the density and pressure force for all the particles
    query.par_iter_mut().for_each(|(mut particle, _)| {
        particle.density = 1.0;
        particle.pressure_force = Vec3::ZERO;
    });


    // Predict the particle positions based on their current velocity
    query.par_iter_mut().for_each(|(mut particle, transform)| {

        particle.position = transform.translation + particle.velocity * time.delta_secs();
    });

    // Get all the particle pairs without repetition
    let mut combinations = query.iter_combinations_mut();

    // Calculate density
    while let Some([(mut particle_a, transform_a), (mut particle_b, transform_b)]) = combinations.fetch_next() {
        // Particle a
        let distance = Vec3::length(particle_b.position - particle_a.position);
        let influence = smoothing_kernel(simulation_parameters.particle_smoothing_radius, distance);
        particle_a.density += influence;

        // Particle b
        let distance = Vec3::length(particle_a.position - particle_b.position);
        let influence = smoothing_kernel(simulation_parameters.particle_smoothing_radius, distance);
        particle_b.density += influence;
    }

    let mut combinations = query.iter_combinations_mut();

    // Calculate pressure force
    while let Some([(mut particle_a, transform_a), (mut particle_b, transform_b)]) = combinations.fetch_next() {

        let shared_pressure_force = (convert_density_to_pressure(particle_b.density, &simulation_parameters) + convert_density_to_pressure(particle_a.density, &simulation_parameters)) / 2.0;

        // Particle a
        let distance = Vec3::length(particle_b.position - particle_a.position);
        let direction = (particle_b.position - particle_a.position) / distance;
        let slope = smoothing_kernel_derivative(simulation_parameters.particle_smoothing_radius, distance);
        particle_a.pressure_force += shared_pressure_force * direction * slope / particle_b.density;
        
        // Particle b
        let distance = Vec3::length(particle_a.position - particle_b.position);
        let direction = (particle_a.position - particle_b.position) / distance;
        let slope = smoothing_kernel_derivative(simulation_parameters.particle_smoothing_radius, distance);
        particle_b.pressure_force += shared_pressure_force * direction * slope / particle_a.density;
    }

    // Apply gravity to all the particles
    /* query.par_iter_mut().for_each(|(mut particle, transform)| {

        particle.velocity += Vec3::NEG_Y * simulation_parameters.gravity * time.delta_secs();
    }); */

    // Apply pressure force to all the particles
    query.par_iter_mut().for_each(|(mut particle, mut transform)| {

        let pressure_acceleration = particle.pressure_force / particle.density;
        particle.velocity += pressure_acceleration * time.delta_secs();
        transform.translation += particle.velocity * time.delta_secs();
    });

}
