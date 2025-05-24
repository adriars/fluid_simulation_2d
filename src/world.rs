use bevy::prelude::*;

pub struct WorldSetup;

impl Plugin for WorldSetup {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_world);
    }
}

// Sets up the world basic configuration
fn setup_world(commands: Commands) {
    spawn_light(commands);
}

// Spawns the main directional light of the world
fn spawn_light(mut commands: Commands) {
    // Light
    commands.spawn(DirectionalLight::default());
}
