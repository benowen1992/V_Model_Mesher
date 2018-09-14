function [ particlePositionX, particlePositionY, num_particles, particle_diameter ] = beam_particle_positions(resolution, height, width, start_x, start_y)
% Beam_particle_positions 
%   Returns particle x and y positions for beam geometry based upon height
%   and width of beam and required resolution

    % Calculate particle diameter for triangle arrangement
    particle_diameter = width / (1 + (resolution - 1) * (0.75)^0.5);
    
    % Calculate number of particle in each direction and in total
    num_y           = resolution;
    num_x           = ceil(height/particle_diameter);
    num_particles   = num_x * num_y;
    
    % Set up particle position vectors 
    particlePositionX = zeros(num_particles,1);                         % particle x positions                       
    particlePositionY = zeros(num_particles,1);                         % particle y positions
     
    % find position of particles
%     for i = 1:num_particles
% 
%         % find row and column value for each particle
%         col_particle = ceil(i / num_x);                                             % calculates particle position in its row (i.e. 1,2,3...)
%         row_particle = i - ((col_particle-1) * num_x);                              % calculates particle position in its column (i.e. 1,2,3...)
% 
%         % find particle y position depending on if in even or odd column
%         if mod(col_particle,2) == 0
%             particlePositionY(i,1) = (row_particle - 1) * particle_diameter + particle_diameter*cos(pi/3);        
%         else 
%             particlePositionY(i,1) = (row_particle - 1) * particle_diameter ;                           
%         end
% 
%         % find particle x position
%         particlePositionX(i,1) = start_x + (col_particle - 1) * particle_diameter * sin(pi/3) + particle_diameter/2 ;                         
%     end
    
    
    % alternative indexing (by row rather than column
    
     % find position of particles
    for i = 1:num_particles

        % find row and column value for each particle
        row_particle = ceil(i/resolution)
        col_particle = i - ((row_particle - 1) * resolution)
        
       % col_particle = ceil(i / num_x);                                             % calculates particle position in its row (i.e. 1,2,3...)
       % row_particle = i - ((col_particle-1) * num_x);                              % calculates particle position in its column (i.e. 1,2,3...)

        % find particle x position 
        % find particle y position depending on if in even or odd column
        if mod(col_particle,2) == 0
            particlePositionY(i,1) = (row_particle - 1) * particle_diameter + particle_diameter*cos(pi/3);        
        else 
            particlePositionY(i,1) = (row_particle - 1) * particle_diameter ;                           
        end

        % find particle x position
        particlePositionX(i,1) = start_x + (col_particle - 1) * particle_diameter * sin(pi/3) + particle_diameter/2 ;                         
    end
    
        
    
    
end

