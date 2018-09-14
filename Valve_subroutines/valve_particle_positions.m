function [ particlePositionX, particlePositionY, num_particles, base_num_particles, leaflet_num_particles, base_particle_diameter, leaflet_particle_diameter, particle_diameter ] = valve_particle_positions( angle, base_width, base_height, leaflet_width, leaflet_height, resolution, start_x, start_y )
%valve_particle_positions 
% find the postions of all particles within idealised valve geometry of
% type Nasar 2016.  Finds particles in base, then leaflet and transforms
% them so they are connected.  Overlapping particles are then deleted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Base particle set-up %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_particle_diameter = base_width / resolution;

base_num_x = resolution;
base_num_y = round(base_height / base_particle_diameter);
base_num_particles = base_num_x * base_num_y;

particlePositionX = zeros(base_num_particles,1);
particlePositionY = zeros(base_num_particles,1);


for i = 1:base_num_particles

    row_particle = ceil(i / base_num_x);                                             % calculates particle position in its row (i.e. 1,2,3...)
    col_particle = i - ((row_particle-1) * base_num_x);                              % calculates particle position in its column (i.e. 1,2,3...)

    particlePositionX(i,1) = start_x + col_particle * base_particle_diameter - base_particle_diameter/2;          % postion of each particle in x direction
    particlePositionY(i,1) = start_y + (row_particle - 1) * base_particle_diameter;                               % postion of each particle in y direction

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Leaflet particle set-up %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Triangular arrangement 

leaflet_particle_diameter = leaflet_width / (1 + (resolution - 1) * (0.75)^0.5);

leaflet_num_x = resolution;
leaflet_num_y = floor(leaflet_height / leaflet_particle_diameter);

leaflet_num_particles = leaflet_num_x * leaflet_num_y;

num_particles = base_num_particles + leaflet_num_particles;


%%% Get particle arrangement in datum positions

 for i = (1 + base_num_particles):(base_num_particles + leaflet_num_particles)

    row_particle = ceil((i - base_num_particles) / leaflet_num_x);                                            % calculates particle position in its row (i.e. 1,2,3...)
    col_particle = (i - base_num_particles) - ((row_particle-1) * leaflet_num_x);                            % calculates particle position in its column (i.e. 1,2,3...)

    if mod(col_particle,2) == 0

        particlePositionY(i,1) = (row_particle - 1) * leaflet_particle_diameter + leaflet_particle_diameter*cos(pi/3);         % postion of each particle in x direction

    else 

         particlePositionY(i,1) = (row_particle - 1) * leaflet_particle_diameter ;                             % postion of each particle in x direction

    end

     particlePositionX(i,1) = (col_particle - 1) * leaflet_particle_diameter * cos(pi/6) + leaflet_particle_diameter/2;   % replace height_beam/2 with more accurate value of particle height from abouzied thesis                            % postion of each particle in y direction

 end

%%% Transform particle arrangement so connected to base

% Translate so particle 1 is located at origin

dx = particlePositionX(base_num_particles + 1,1);
dy = particlePositionY(base_num_particles + 1,1);

for i = (1 + base_num_particles):num_particles
    
    particlePositionX(i,1) = particlePositionX(i,1) - dx;
    particlePositionY(i,1) = particlePositionY(i,1) - dy;
    
end

% Rotate to angle

for i = (1 + base_num_particles):num_particles
    
    temp_x = particlePositionX(i,1) * cos(angle) - particlePositionY(i,1) * sin(angle);
    temp_y = particlePositionX(i,1) * sin(angle) + particlePositionY(i,1) * cos(angle);
    
    particlePositionX(i,1) = temp_x;
    particlePositionY(i,1) = temp_y;
    
end

% Translate to base

base_datum_x = start_x + base_particle_diameter/2;
base_datum_y = base_height;

for i = (1 + base_num_particles):num_particles
    
    particlePositionX(i,1) = particlePositionX(i,1) + base_datum_x;
    particlePositionY(i,1) = particlePositionY(i,1) + base_datum_y;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Connection clean up %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

particles_delete = zeros(1,1);

for i = (1 + base_num_particles):num_particles
    
    if particlePositionY(i,1) < base_height
        
        particles_delete = [particles_delete, i];
    end
    
end

particles_delete(1) = [];
particles_delete_reverse = fliplr(particles_delete);

for i = 1:length(particles_delete)
    
    particlePositionX(particles_delete_reverse(i)) = [];   
    particlePositionY(particles_delete_reverse(i)) = [];    

end


num_particles           = length(particlePositionX);
leaflet_num_particles   = num_particles - base_num_particles;
particle_diameter       = leaflet_particle_diameter;


end

