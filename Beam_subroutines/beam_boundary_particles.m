function [ particle_boundary, num_boundary_particles, x_desired, y_desired  ] = beam_boundary_particles( height, width, start_x, start_y, particle_diameter, particlePositionX, particlePositionY, num_particles )
%beam_boundary_particles 
%   Identifies the particles at the boundary of the beam by searching for
%   particles within 1 particle diameter of the desired geometry


% Corner coordinates of desired gemoetry
x_desired = [start_x start_x start_x + width start_x + width];
y_desired = [start_y start_y + height start_y + height start_y];
    
%%% Discretise into small points for nearest neighbour search
points = round(width/particle_diameter * 50);

geom_x = [linspace(x_desired(1),x_desired(2),points) linspace(x_desired(2),x_desired(3),points) linspace(x_desired(3),x_desired(4),points)];
geom_y = [linspace(y_desired(1),y_desired(2),points) linspace(y_desired(2),y_desired(3),points) linspace(y_desired(3),y_desired(4),points)];

particle_boundary = zeros(1,1);

%%% Loop through desired geometry points to find boundary particles

for i = 1:length(geom_x)
    for ii = 1:num_particles
        
        dist_mag = sqrt((geom_x(i)-particlePositionX(ii))^2 + (geom_y(i)-particlePositionY(ii))^2) ;

  
        
        if abs(dist_mag) < 0.7 * particle_diameter 

            particle_boundary = [particle_boundary ii];

        end
    end
end

particle_boundary(1)    = [];
particle_boundary       = unique(particle_boundary);
num_boundary_particles  = length(particle_boundary);

end

