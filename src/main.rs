use std::collections::HashMap;
use raylib::prelude::*;
use rand::prelude::*;

type ParticleKindId = usize;
struct ParticleKind {
    color: Color,
    radius: f32,
}

type ParticleId = usize;
#[derive(Clone)]
struct Particle {
    pos: Vector2,
    vel: Vector2,
    kind_id: ParticleKindId,
}

#[derive(Default)]
struct Simulation {
    size: (f32, f32),

    force_coefficients: HashMap<(ParticleKindId, ParticleKindId), f32>,

    next_particle_id: ParticleId,
    particles: HashMap<ParticleId, Particle>,

    next_particle_kind_id: ParticleKindId,
    particle_kinds: HashMap<ParticleKindId, ParticleKind>,
}

impl Simulation {
    pub fn new() -> Self { Self::default() }

    pub fn add_random_kind(&mut self) -> ParticleKindId {
        let color = Color::color_from_hsv(rand::thread_rng().gen_range(0.0..360.0), 0.9, 0.9);
        let radius = 15.0;
        let kind = ParticleKind { color, radius };
        let kind_id = self.next_particle_kind_id;

        for &other_kind_id in self.particle_kinds.keys() {
            let coefficient = rand::thread_rng().gen_range(-1.0..1.0);
            self.force_coefficients.insert((kind_id, other_kind_id), coefficient);
            let coefficient = rand::thread_rng().gen_range(-1.0..1.0);
            self.force_coefficients.insert((other_kind_id, kind_id), coefficient);
        }
        let coefficient = rand::thread_rng().gen_range(-1.0..1.0);
        self.force_coefficients.insert((kind_id, kind_id), coefficient);
        self.next_particle_kind_id += 1;
        self.particle_kinds.insert(kind_id, kind);
        kind_id
    }

    pub fn add_random_particle(&mut self) -> ParticleId {
        let kind_id = *self.particle_kinds.keys().choose(&mut rand::thread_rng()).unwrap();
        let pos = Vector2::new(
            rand::thread_rng().gen_range(-self.size.0..self.size.0),
            rand::thread_rng().gen_range(-self.size.1..self.size.1),
        );
        let vel = Vector2::zero();
        let particle = Particle { pos, vel, kind_id };
        let particle_id = self.next_particle_id;

        self.next_particle_id += 1;
        self.particles.insert(particle_id, particle);
        particle_id
    }

    pub fn draw<Draw: RaylibDraw>(&self, d: &mut Draw) {
        d.clear_background(Color::BLACK);
        for Particle { pos, kind_id, .. } in self.particles.values() {
            let ParticleKind { color, radius } = self.particle_kinds[kind_id];
            d.draw_circle_v(pos, radius, color);
        }
    }

    pub fn step(&mut self, delta_time: f32) {
        let old_particles = self.particles.clone();
        for (&id, particle) in &mut self.particles {
            let particle_kind = &self.particle_kinds[&particle.kind_id];

            particle.pos += old_particles[&id].vel * delta_time;
            let (max_x, max_y) = (self.size.0 - particle_kind.radius, self.size.1 - particle_kind.radius);
            if particle.pos.x.abs() > max_x {
                particle.pos.x = max_x * particle.pos.x.signum();
                particle.vel.x *= -1.0;
            }
            if particle.pos.y.abs() > max_y {
                particle.pos.y = max_y * particle.pos.y.signum();
                particle.vel.y *= -1.0;
            }
            particle.vel *= 0.99;
            for (&other_id, other_particle) in &old_particles {
                let other_particle_kind = &self.particle_kinds[&other_particle.kind_id];

                if id != other_id {
                    let dir = other_particle.pos - particle.pos;
                    let dist = dir.length();
                    let dir = dir / dist;
                    if dist > self.particle_kinds[&particle.kind_id].radius + self.particle_kinds[&other_particle.kind_id].radius {
                        let coefficient = self.force_coefficients[&(particle.kind_id, other_particle.kind_id)];
                        particle.vel += dir / dist / dist * coefficient * 2000000.0 * delta_time;
                    } else {
                        particle.vel -= dir * 20.0 * delta_time * 10.0;
                    }
                }
            }
        }
    }
}


fn main() {
    let mut sim = Simulation::new();
    sim.size = (640.0, 480.0);
    for _ in 0..5 { sim.add_random_kind(); };
    for _ in 0..150 { sim.add_random_particle(); };
    let (mut rl, thread) = raylib::init()
        .size(640, 480)
        .title("Hello, World")
        .build();

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::WHITE);
        let mut d = d.begin_mode2D(Camera2D {
            offset: Vector2::from(sim.size) * 0.5,
            target: Vector2::zero(),
            rotation: 0.0,
            zoom: 0.5,
        });
        sim.draw(&mut d);
        sim.step(d.get_frame_time());
    }
}
