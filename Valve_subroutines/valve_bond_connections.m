function  [ b2p_ID, num_bonds ] = valve_bond_connections( num_particles, particlePositionX, particlePositionY, base_particle_diameter, leaflet_particle_diameter, base_height);
% valve_bond_connections - identifies bonds which connect particles
% Loops through particles doing nearest neighbour search. If distance
% between particle is around size of particle diameter then a bond is
% formed between them. Bonds are then sorted so only one bond exists
% between any two particles.
% Loops through base particles and leaflet particles seperately since they
% have different bond lengths.  May require interation to create bonds
% at interface between base and leaflet sections

%%% Find boundary bonds

% create bonds between all particles 

b2p_IDall = zeros(1,2);

 for i = 1:num_particles
    for ii = 1:num_particles

        % Do nearest neighbour search
        rx_dist = particlePositionX(ii,1) - particlePositionX(i,1);
        ry_dist = particlePositionY(ii,1) - particlePositionY(i,1);
        Mag_r_dist = sqrt(rx_dist^2 + ry_dist^2);

        % Find bonds in base
        if particlePositionY(ii,1) < base_height * 0.99 && Mag_r_dist <= 1.08 * base_particle_diameter && Mag_r_dist >= 0.5 * base_particle_diameter
            bond = [i ii];
            b2p_IDall = [b2p_IDall;bond];
        end
            
        % Find leaflets in base
        if particlePositionY(ii,1) > base_height * 0.99 && Mag_r_dist <= 1.1 * leaflet_particle_diameter && Mag_r_dist >= 0.5 * leaflet_particle_diameter
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

