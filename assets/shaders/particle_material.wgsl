#import bevy_pbr::forward_io::VertexOutput

@group(0) @binding(0) var<uniform> cell_key: u32;

@fragment
fn fs_main(mesh: VertexOutput) -> @location(0) vec4<f32> {

    return vec4<f32>(f32(cell_key), f32(cell_key), f32(cell_key), 1.0);
}