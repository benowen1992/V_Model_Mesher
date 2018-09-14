function [b2p_ID_boundary_sequence, num_boundary_bonds] = boundary_bonds(b2p_ID, particle_boundary, num_boundary_particles, particlePositionX, particlePositionY, start_x, start_y)
% boundary bonds 
%   Loops through the boundary particles to find the boundary bonds in
%   sequence to place IB markers along them with spacing equal to the grid
%   resolution.


%%% find bonds between boundary particles

% find boundary bonds
b2p_ID_boundary = zeros(1,2);

for i = 1:length(particle_boundary)
    for ii = 1:length(particle_boundary)
        for iii = 1:length(b2p_ID)
            if b2p_ID(iii,1) == particle_boundary(i) && b2p_ID(iii,2) == particle_boundary(ii)
                 b2p_ID_boundary = [b2p_ID_boundary; b2p_ID(iii,:)];
            end
        end
    end
end

b2p_ID_boundary(1,:) = [];
b2p_ID_boundary;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% sequence boundary bonds %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find first bond relative to start position

min_1 = 1000;
min_2 = 10000;

start_particle_i = 1000;
start_particle_j = 10000;


    for i = 1:num_boundary_particles

            % Do nearest neighbour search
            rx_dist = particlePositionX(i) - start_x;
            ry_dist = particlePositionY(i) - start_y;
            Mag_r_dist = sqrt(rx_dist^2 + ry_dist^2);

            if Mag_r_dist < min_1

                min_2 = min_1;
                min_1 = Mag_r_dist;

                start_particle_j = start_particle_i;
                start_particle_i = i;

            end

            if Mag_r_dist < min_2 && Mag_r_dist > min_1

                min_2 = Mag_r_dist;

                start_particle_j = i;

            end
    end
        
start_bond = [start_particle_i start_particle_j];

b2p_ID_boundary_sequence = start_bond;
previous_bond = start_bond;
connecting_particle = start_bond(2);
last_connecting_particle = start_bond(1);

i = 1;
ii = 1;

    while length(b2p_ID_boundary_sequence(:,1)) < length(b2p_ID_boundary(:,1))

        if b2p_ID_boundary(i,1) == connecting_particle && b2p_ID_boundary(i,2) ~= last_connecting_particle

            b2p_ID_boundary_sequence = [b2p_ID_boundary_sequence; b2p_ID_boundary(i,:)];

            last_connecting_particle = connecting_particle;
            connecting_particle      = b2p_ID_boundary(i,2);

        end

         if b2p_ID_boundary(i,2) == connecting_particle && b2p_ID_boundary(i,1) ~= last_connecting_particle

            b2p_ID_boundary_sequence = [b2p_ID_boundary_sequence; b2p_ID_boundary(i,2) b2p_ID_boundary(i,1) ];
            last_connecting_particle = connecting_particle;
            connecting_particle      = b2p_ID_boundary(i,1);

         end

         i = i + 1;

         if i > length(b2p_ID_boundary(:,1))
             i = 1;
             ii = ii + 1;
         end

         if ii > 500
             break
         end
    end
    
    num_boundary_bonds = length(b2p_ID_boundary_sequence);
end

