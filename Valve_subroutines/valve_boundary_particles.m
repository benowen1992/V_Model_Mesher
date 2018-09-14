function [ particle_boundary, num_boundary_particles, x_desired, y_desired ] = valve_boundary_particles( base_height, base_width, angle, leaflet_height, leaflet_width, start_x, start_y, base_particle_diameter, leaflet_particle_diameter, particlePositionX, particlePositionY, num_particles)
%valve_boundary_particles 
% Finds particles on the boundary of the desired shape by finding particles
% within one particle diameter of the desired geoemetry around the
% parameter.

% Corner coordinates of desired gemoetry

x0 = start_x;
x1 = x0;
x2 = start_x + leaflet_height * cos(pi/2 - abs(angle));
x3 = x2 + leaflet_width * cos(abs(angle));
x4 = start_x + base_width;
x5 = x4;

y0 = start_y;
y1 = start_y + base_height;
y2 = y1 + leaflet_height * sin(pi/2 - abs(angle));
y3 = y2 - leaflet_width * sin(abs(angle));
y4 = y1;
y5 = y0;


x_desired = [x0 x1 x2 x3 x4 x5];
y_desired = [y0 y1 y2 y3 y4 y5];


%%% Discretise into small points for nearest neighbour search
points = 200;

geom_x = [linspace(x0,x1,points) linspace(x1,x2,points) linspace(x2,x3,points) linspace(x3,x4,points) linspace(x4,x5,points)];
geom_y = [linspace(y0,y1,points) linspace(y1,y2,points) linspace(y2,y3,points) linspace(y3,y4,points) linspace(y4,y5,points)];

particle_boundary = zeros(1,1);

%%% Loop through desired geometry points to find boundary particles

for i = 1:length(geom_x)
    for ii = 1:length(particlePositionX(:,1))
        
        dist_mag = sqrt((geom_x(i)-particlePositionX(ii,1))^2 + (geom_y(i)-particlePositionY(ii,1))^2) ;

        if particlePositionY(ii,1) < y1 * 1.01 && abs(dist_mag) < 0.7 * base_particle_diameter 

            particle_boundary = [particle_boundary ii];

        end
        
        if particlePositionY(ii,1) > y1 && abs(dist_mag) < 1.0 * leaflet_particle_diameter 

            particle_boundary = [particle_boundary ii];

        end
    end
end

particle_boundary(1) = [];
particle_boundary = unique(particle_boundary);

num_boundary_particles = length(particle_boundary);





end

