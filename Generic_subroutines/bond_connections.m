function [ b2p_ID, num_bonds ] = bond_connections( num_particles, particlePositionX, particlePositionY, particle_diameter)
%bond_connections - identifies bonds which connect particles
%   Loops through particles doing nearest neighbour search. If distance
%   between particle is around size of particle diameter then a bond is
%   formed between them. Bonds are then sorted so only one bond exists
%   between any two particles.

% create bonds between all particles 

b2p_IDall = zeros(1,2);

 for i = 1:num_particles
    for ii = 1:num_particles

        % Do nearest neighbour search
        rx_dist = particlePositionX(ii) - particlePositionX(i);
        ry_dist = particlePositionY(ii) - particlePositionY(i);
        Mag_r_dist = sqrt(rx_dist^2 + ry_dist^2);
            
        % Check if particles are close enough for bond between them
        if  Mag_r_dist <= 1.1 * particle_diameter && Mag_r_dist >= 0.9 * particle_diameter
            bond = [i ii];
            b2p_IDall = [b2p_IDall;bond];
        end
    end
 end
 
 b2p_IDall;
             
% Sort bond to particle arrays     
b2p_IDall(1,:)=[];
b2p_IDall;
b2p_IDsorted = b2p_IDall;
      
for i = 1:length(b2p_IDall)
    if b2p_IDall(i,2) < b2p_IDall(i,1)
        b2p_IDsorted(i,2) = b2p_IDall(i,1);
        b2p_IDsorted(i,1) = b2p_IDall(i,2);
    end
end
    
b2p_ID      = unique(b2p_IDsorted,'rows');   
num_bonds   = length(b2p_ID);


end

