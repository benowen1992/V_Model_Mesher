function [] = valve_figure_output(particlePositionX, particlePositionY, num_particles, base_num_particles, leaflet_num_particles, base_particle_diameter, leaflet_particle_diameter, num_boundary_particles, particle_boundary, start_x, start_y, leaflet_height, x_desired, y_desired, num_bonds, b2p_ID, b2p_ID_boundary_sequence, num_boundary_bonds, num_clamped_particles, clamped_particles)
% valve figure output 
% Outputs figure to check set up for:
% 1) Particle positions and desired geometry
% 2) 
% 3)

figure(1)
rad_global = zeros(num_particles,1);
rad_global(1:base_num_particles) = base_particle_diameter/2;
rad_global(base_num_particles+1 : num_particles) = leaflet_particle_diameter/2;
r = [particlePositionX particlePositionY];
viscircles(r,rad_global,'LineWidth',2);
axis equal
ylim([start_y leaflet_height])
hold on

%%% desired geometry
plot(x_desired, y_desired, 'k')


% bonds connections
figure(2)
r = [particlePositionX particlePositionY];
viscircles(r,rad_global,'LineWidth',2);
axis equal
ylim([start_y leaflet_height])
hold on 

for i = 1:num_bonds
    x_bond = [particlePositionX(b2p_ID(i,1)) particlePositionX(b2p_ID(i,2))];
    y_bond = [particlePositionY(b2p_ID(i,1)) particlePositionY(b2p_ID(i,2))];

    plot(x_bond,y_bond,'k','LineWidth',1);
end



% boundary particles

boundary_particlePositionX = zeros(num_boundary_particles,1);
boundary_particlePositionY = zeros(num_boundary_particles,1);
rad_boundary = zeros(num_boundary_particles,1);


for i = 1:num_boundary_particles

    boundary_particlePositionX(i) = particlePositionX(particle_boundary(i));
    boundary_particlePositionY(i) = particlePositionY(particle_boundary(i));
    rad_boundary(i)               = rad_global(particle_boundary(i));

end

r = [boundary_particlePositionX boundary_particlePositionY];



figure(3)
viscircles(r,rad_boundary,'LineWidth',2);
axis equal
ylim([start_y leaflet_height])
hold on


figure(4)
for i = 1:length(b2p_ID_boundary_sequence(:,1))
    
    x_bond = [particlePositionX(b2p_ID_boundary_sequence(i,1),1) particlePositionX(b2p_ID_boundary_sequence(i,2),1)];
    y_bond = [particlePositionY(b2p_ID_boundary_sequence(i,1),1) particlePositionY(b2p_ID_boundary_sequence(i,2),1)];
    
    plot(x_bond,y_bond,'k');
    hold on
end
axis equal
ylim([start_y leaflet_height])

figure(5)
rad_global = zeros(num_particles,1);
rad_global(1:base_num_particles) = base_particle_diameter/2;
rad_global(base_num_particles+1 : num_particles) = leaflet_particle_diameter/2;
r = [particlePositionX particlePositionY];
viscircles(r,rad_global,'LineWidth',2);
axis equal
ylim([start_y leaflet_height])
hold on

particleClampedPositionX = zeros(num_clamped_particles,1);
particleClampedPositionY = zeros(num_clamped_particles,1);
rad_clamped = zeros(num_clamped_particles,1);


for i = 1:num_clamped_particles
   
    particleClampedPositionX(i) = particlePositionX(clamped_particles(i));
    particleClampedPositionY(i) = particlePositionY(clamped_particles(i));
    rad_clamped(i)              = rad_global(clamped_particles(i));
    
end

r_clamped = [particleClampedPositionX particleClampedPositionY];
viscircles(r_clamped,rad_clamped,'EdgeColor','k','LineWidth',2);



