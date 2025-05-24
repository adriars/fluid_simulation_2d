@group(2) @binding(0)
var<uniform> velocity: f32;

@fragment
fn fragment() -> @location(0) vec4<f32> {
  return vec4<f32>(velocity, velocity, velocity, 1.0);
}