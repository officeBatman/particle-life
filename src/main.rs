use std::collections::HashMap;
use raylib::prelude::*;
use rand::prelude::*;

fn random_on_circle<R: RngCore>(rng: &mut R) -> Vector2 {
    let angle: f32 = rng.gen_range(0.0..=360f32.to_radians());
    Vector2::new(angle.cos(), angle.sin())
}

fn random_in_circle<R: RngCore>(rng: &mut R) -> Vector2 {
    let point = Vector2::new(rng.gen_range(-1.0..=1.0), rng.gen_range(-1.0..=1.0));
    if point.length_sqr() <= 1.0 { point } else { random_in_circle(rng) }
}

type ParticleKindId = usize;
struct ParticleKind {
    color: Color,
    radius: f32,
    force_max_distance: f32,
}

type ParticleId = usize;
#[derive(Clone)]
struct Particle {
    pos: Vector2,
    vel: Vector2,
    kind_id: ParticleKindId,
}

struct SimulationConfig {
    size: (f32, f32),
    rng: Box<dyn RngCore>,
    particle_coloring: fn (&mut dyn RngCore) -> Color,
    particle_radius: fn (&mut dyn RngCore) -> f32,
    generate_force_max_distance: fn (&mut dyn RngCore) -> f32,
    force_constant: f32,
    energy: f32,
    max_vel: f32,
    flacuation_constant: f32,
    collision_vel_damping: f32,
}

struct Simulation {
    cfg: SimulationConfig,

    force_coefficients: HashMap<(ParticleKindId, ParticleKindId), f32>,

    next_particle_id: ParticleId,
    particles: HashMap<ParticleId, Particle>,

    next_particle_kind_id: ParticleKindId,
    particle_kinds: HashMap<ParticleKindId, ParticleKind>,
}

impl Simulation {
    pub fn new(cfg: SimulationConfig) -> Self {
        Self {
            cfg,
            force_coefficients: Default::default(),
            next_particle_id: 0,
            particles: Default::default(),
            next_particle_kind_id: 0,
            particle_kinds: Default::default(),
        }
    }

    pub fn size(&self) -> (f32, f32) { self.cfg.size }
    pub fn size_mut(&mut self) -> &mut (f32, f32) { &mut self.cfg.size }

    pub fn width(&self) -> f32 { self.cfg.size.0 }
    pub fn height(&self) -> f32 { self.cfg.size.1 }
    
    pub fn rng(&self) -> &dyn RngCore { &self.cfg.rng }
    pub fn rng_mut(&mut self) -> &mut dyn RngCore { &mut self.cfg.rng }

    pub fn add_random_kind(&mut self) -> ParticleKindId {
        let color = (self.cfg.particle_coloring)(self.rng_mut());
        let radius = (self.cfg.particle_radius)(self.rng_mut());
        let force_max_distance = (self.cfg.generate_force_max_distance)(self.rng_mut());
        let kind = ParticleKind { color, radius, force_max_distance };
        let kind_id = self.next_particle_kind_id;

        for &other_kind_id in self.particle_kinds.keys() {
            let coefficient = self.cfg.rng.gen_range(-1.0..1.0);
            self.force_coefficients.insert((kind_id, other_kind_id), coefficient);
            let coefficient = self.cfg.rng.gen_range(-1.0..1.0);
            self.force_coefficients.insert((other_kind_id, kind_id), coefficient);
        }
        let coefficient = self.rng_mut().gen_range(-1.0..1.0);
        self.force_coefficients.insert((kind_id, kind_id), coefficient);
        self.next_particle_kind_id += 1;
        self.particle_kinds.insert(kind_id, kind);
        kind_id
    }

    pub fn add_random_particle(&mut self) -> ParticleId {
        let kind_id = self.particle_kinds.keys().choose(&mut self.cfg.rng).cloned().unwrap();
        let pos = Vector2::new(
            self.cfg.rng.gen_range(-0.5*self.width()..0.5*self.width()),
            self.cfg.rng.gen_range(-0.5*self.height()..0.5*self.height()),
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
            let ParticleKind { color, radius, .. } = self.particle_kinds[kind_id];
            d.draw_circle_v(pos, radius, color);
        }
    }

    fn border_collision(particle: &mut Particle, particle_kind: &ParticleKind, sim_size: (f32, f32)) {
        let max_x = (sim_size.0 - particle_kind.radius) * 0.5;
        let max_y = (sim_size.1 - particle_kind.radius) * 0.5;

        if particle.pos.x.abs() > max_x {
            particle.pos.x = max_x * particle.pos.x.signum();
            particle.vel.x *= -1.0;
        }
        if particle.pos.y.abs() > max_y {
            particle.pos.y = max_y * particle.pos.y.signum();
            particle.vel.y *= -1.0;
        }
    }

    fn distribute_energy(particle_vel: &mut Vector2, total_energy: f32, target_energy: f32) {
        if total_energy == target_energy || total_energy == 0.0 { return; }

        *particle_vel *= (target_energy / total_energy).sqrt();
    }

    fn force_formula(k: f32, dist: f32, min_dist: f32, max_dist: f32) -> f32 {
        k * (max_dist - dist) / (max_dist - min_dist)
    }

    pub fn step(&mut self, delta_time: f32) {
        let old_particles = self.particles.clone();
        for (&id, particle) in &mut self.particles {
            let particle_kind = &self.particle_kinds[&particle.kind_id];

            //  Move the particle
            particle.pos += old_particles[&id].vel * delta_time;
            //  Collision
            Self::border_collision(particle, particle_kind, self.cfg.size);
            //  Flacuation
            particle.vel += random_in_circle(&mut self.cfg.rng) * self.cfg.flacuation_constant * delta_time;

            for (&other_id, other_particle) in &old_particles {
                if id != other_id {
                    let other_particle_kind = &self.particle_kinds[&other_particle.kind_id];

                    let dir = other_particle.pos - particle.pos;
                    let square_dist = dir.length_sqr();
                    let min_dist = particle_kind.radius + other_particle_kind.radius;
                    let max_dist = particle_kind.force_max_distance; 
                    if square_dist > max_dist * max_dist || square_dist == 0. {
                        //  nothing
                    } else {
                        let dist = dir.length();
                        let dirn = dir / dist;
                        if dist > min_dist {
                            let coefficient = self.force_coefficients[&(particle.kind_id, other_particle.kind_id)];
                            let force = Self::force_formula(self.cfg.force_constant, dist, particle_kind.radius, particle_kind.force_max_distance);
                            particle.vel += dirn * force * coefficient * delta_time;
                        } else {
                            //  move the rigidbody half way out - the other half is moved in
                            //  another iteration of the loop
                            particle.pos -= dirn * (min_dist - dist) / 2.;
                            particle.vel -= dirn * particle.vel.dot(dirn) * self.cfg.collision_vel_damping;
                        }
                    }
                    if particle.vel.length_sqr() > self.cfg.max_vel * self.cfg.max_vel {
                        particle.vel = particle.vel.normalized() * self.cfg.max_vel;
                    }
                }
            }
        }

        let total_energy = self.particles.values().map(|p| p.vel.length_sqr()).sum();
        for particle in self.particles.values_mut() {
            Self::distribute_energy(&mut particle.vel, total_energy, self.cfg.energy);
        }
    }
}


fn main() {
    let (mut rl, thread) = raylib::init()
        .fullscreen()
        .size(0, 0)
        .title("Hello, World")
        .build();

    let mut sim = Simulation::new(SimulationConfig {
        size: (rl.get_screen_width() as f32 / rl.get_screen_height() as f32 * 100., 100.),
        rng: Box::new(rand::thread_rng()),
        particle_coloring: |rng| Color::color_from_hsv(rng.gen_range(0.0..360.0), 0.9, 0.9),
        particle_radius: |rng| rng.gen_range(0.8..=1.2),
        generate_force_max_distance: |rng| rng.gen_range(1.0..60.0),
        force_constant: 50.0,
        energy: 20000.0,
        max_vel: 30.0,
        flacuation_constant: 20.0,
        collision_vel_damping: 0.5,
    });
    for _ in 0..10 { sim.add_random_kind(); };
    for _ in 0..300 { sim.add_random_particle(); };

    while !rl.window_should_close() {
        let camera = Camera2D {
            offset: Vector2::new(rl.get_screen_width() as _, rl.get_screen_height() as _) * 0.5,
            target: Vector2::zero(),
            rotation: 0.0,
            zoom: rl.get_screen_width() as f32 / sim.cfg.size.0,
        };
        let mut screen_draw = rl.begin_drawing(&thread);
        let mut world_draw = screen_draw.begin_mode2D(camera);
        sim.draw(&mut world_draw);
        sim.step(world_draw.get_frame_time());
    }
}
