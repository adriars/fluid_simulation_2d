mod camera;
mod world;
mod fluid_simulation;
mod bounding_box;
mod ui_debug;

use bevy::{
    app::App,
    DefaultPlugins,
};
use bounding_box::BoundingBox;
use camera::MainCamera;
use fluid_simulation::FluidSimulation;
use world::WorldSetup;
use ui_debug::UIDebug;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(MainCamera)
        .add_plugins(WorldSetup)
        .add_plugins(FluidSimulation)
        .add_plugins(BoundingBox)
        .add_plugins(UIDebug)
        .run();
}
