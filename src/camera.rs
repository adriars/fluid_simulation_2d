use bevy::{app::{App, Plugin, PostStartup, PreUpdate, Startup, Update}, core_pipeline::core_2d::Camera2d, ecs::{component::Component, query::With, system::{Commands, Query, Res, Single}}, input::{keyboard::KeyCode, mouse::AccumulatedMouseScroll, ButtonInput}, render::camera::{Camera, OrthographicProjection, Projection}, ui::widget::Text, utils::default};

pub struct MainCamera;

impl Plugin for MainCamera {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, spawn_camera)
        .add_systems(Update, zoom);
    }
}

// Spawns the main camera of the scene
fn spawn_camera(
    mut commands: Commands
) {
    commands.spawn((
        Camera2d,
        Projection::from(OrthographicProjection {
            scale: 0.05,
            ..OrthographicProjection::default_2d()
        })
    )
    );
}

fn zoom(
    camera: Single<&mut Projection, With<Camera2d>>,
    mouse_wheel_input: Res<AccumulatedMouseScroll>,
    keyboard_input: Res<ButtonInput<KeyCode>>,
) {

    if keyboard_input.pressed(KeyCode::KeyZ) && keyboard_input.pressed(KeyCode::ShiftLeft) {
    
    /*    // Usually, you won't need to handle both types of projection,
    // but doing so makes for a more complete example.
    match *camera.into_inner() {
        Projection::Orthographic(ref mut orthographic) => {
            // We want scrolling up to zoom in, decreasing the scale, so we negate the delta.
            let delta_zoom = -mouse_wheel_input.delta.y * 0.1;
            // When changing scales, logarithmic changes are more intuitive.
            // To get this effect, we add 1 to the delta, so that a delta of 0
            // results in no multiplicative effect, positive values result in a multiplicative increase,
            // and negative values result in multiplicative decreases.
            let multiplicative_zoom = 1. + delta_zoom;

            orthographic.scale = (orthographic.scale * multiplicative_zoom);
        }
        _ => (),
    } */

    }

    if keyboard_input.pressed(KeyCode::KeyZ) {

        // Usually, you won't need to handle both types of projection,
    // but doing so makes for a more complete example.
    match *camera.into_inner() {
        Projection::Orthographic(ref mut orthographic) => {
            // We want scrolling up to zoom in, decreasing the scale, so we negate the delta.
            let delta_zoom = -mouse_wheel_input.delta.y * 0.1;
            // When changing scales, logarithmic changes are more intuitive.
            // To get this effect, we add 1 to the delta, so that a delta of 0
            // results in no multiplicative effect, positive values result in a multiplicative increase,
            // and negative values result in multiplicative decreases.
            let multiplicative_zoom = 1. + delta_zoom;

            orthographic.scale = (orthographic.scale * multiplicative_zoom);
        }
        _ => (),
    }

    }
}