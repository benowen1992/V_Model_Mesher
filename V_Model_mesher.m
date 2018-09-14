%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V-Model Mesher %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% created  27/3/18 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% by Ben Owen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FSI Simulation %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FSI = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Output files/figures %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_files    = 1;
output_figures  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Geometry type flags %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam        = 1;         % a rectangle of uniform length and width
valve       = 0;         % a rigid base with flexible beam attached to top at given angle
hron_turek  = 0;         % a rigid circle with flexible beam attached aft of circle 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Boundary Condition  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clamped_at_wall = 1;           % all of geometry is flexible apart from particles in contact with wall
clamped_at_base = 0;           % some of geometry is rigid between wall and flexible part


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Add subpaths %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Generic_subroutines/')

if beam == 1
    addpath('Beam_subroutines/')
end

if valve == 1
    addpath('Valve_subroutines/')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Particle model parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resolution = 5;         % this usually refers to the number of particles across width           


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Geometry dimensions %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if beam == 1

    height  = 1.0;        % y direction
    width   = 0.06;      % x direction
    start_x = 1.0;        % bottom left corner x coordinate
    start_y = 0;        % bottom left corner y coordinate

end

if valve == 1
    
    angle = -pi/3;                      % angle between leaflet and base

    % leaflet dimensions
    leaflet_height = 0.00844;
    leaflet_width  = 0.0003;

    % base dimensions
    base_height = 0.00141;
    base_width  = leaflet_width / sin(pi/2 - abs(angle));
    
    start_x = 0.0225;        % bottom left corner x coordinate
    start_y = 0;        % bottom left corner y coordinate

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Boundary conditions %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if clamped_at_wall == 1
       
    clamped_postion = start_y;              % only particles connected to wall are clam[ed
    
end


if clamped_at_base == 1
       
    clamped_postion = base_height;          % less than this height is clamped
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Particle arrangement setup %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if beam == 1
    
    [particlePositionX, particlePositionY, num_particles, particle_diameter] = beam_particle_positions(resolution, height, width, start_x, start_y);

end

if valve == 1
   
    [particlePositionX, particlePositionY, num_particles, base_num_particles, leaflet_num_particles, base_particle_diameter, leaflet_particle_diameter, particle_diameter ] = valve_particle_positions( angle, base_width, base_height, leaflet_width, leaflet_height, resolution, start_x, start_y );

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Bond connections %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if beam == 1

    [ b2p_ID, num_bonds ] = beam_bond_connections( num_particles, particlePositionX, particlePositionY, particle_diameter);

end

if valve == 1

    [ b2p_ID, num_bonds ] = valve_bond_connections( num_particles, particlePositionX, particlePositionY, base_particle_diameter, leaflet_particle_diameter, base_height);

end


if FSI == 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Find Boundary Particles %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if beam == 1

        [particle_boundary, num_boundary_particles, x_desired, y_desired ] = beam_boundary_particles( height, width, start_x, start_y, particle_diameter, particlePositionX, particlePositionY, num_particles );

    end

    if valve == 1

        [particle_boundary, num_boundary_particles, x_desired, y_desired ] = valve_boundary_particles( base_height, base_width, angle, leaflet_height, leaflet_width, start_x, start_y, base_particle_diameter, leaflet_particle_diameter, particlePositionX, particlePositionY, num_particles );

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Find and Sequence Boundary Bonds %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [b2p_ID_boundary_sequence, num_boundary_bonds] = boundary_bonds( b2p_ID, particle_boundary, num_boundary_particles, particlePositionX, particlePositionY, start_x, start_y);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set Boundary conditions on particles %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%%% clamped %%%
%%%%%%%%%%%%%%%

if clamped_at_wall == 1
    
   clamped_position  = start_y + 0.8 * particle_diameter;
   clamped_particles = zeros(1,1);
   
   for i = 1:num_particles
      
       if particlePositionY(i) <= clamped_position
       
           clamped_particles = [clamped_particles i];
           
       end
   end
   
   clamped_particles(1) = []; 
   num_clamped_particles = length(clamped_particles);
   
end


if clamped_at_base == 1
    
   clamped_position  = base_height;
   clamped_particles = zeros(1,1);
   
   for i = 1:num_particles
      
       if particlePositionY(i) <= clamped_position
       
           clamped_particles = [clamped_particles i];
           
       end
   end
   
   clamped_particles(1) = []; 
   num_clamped_particles = length(clamped_particles);
   
end


%particlePositionY = -particlePositionY + 0.012;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Print out to text files %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if output_files == 1
    
%%%% Remove previous files and recreate directory to store files    
rmdir Mesh s

mkdir Mesh
mkdir Mesh/Figures


%%%% File with geometry and mesh description
fileID = fopen('Mesh/Geometry and Mesh Description.txt','w');
fprintf(fileID,'%12s\n','Geometry Type:');

    if valve == 1
        type = 'Valve';
    end
    
    if beam == 1
         type = 'Beam';
    end
    
    if clamped_at_wall == 1
         BCtype = ' clamped at wall';
    end
    
    if clamped_at_base == 1
         BCtype = ' with rigid base';
    end
    
    
    
fprintf(fileID,'%s',type);
fprintf(fileID,'%s\n',BCtype);
fprintf(fileID,'\n');
fprintf(fileID,'%12s\n','Number of Particles:');
fprintf(fileID,'%d',num_particles);
fprintf(fileID,'\n');
fprintf(fileID,'%12s\n','Number of Bonds:');
fprintf(fileID,'%d',num_bonds);
fprintf(fileID,'\n');

if FSI == 1
    fprintf(fileID,'%12s\n','Number of Boundary Particles:');
    fprintf(fileID,'%d',num_boundary_particles);
    fprintf(fileID,'\n');
    fprintf(fileID,'%12s\n','Number of Boundary Bonds:');
    fprintf(fileID,'%d',num_boundary_bonds);
end





    
%%%% File with particle and bond properties and numbers
fileID = fopen('Mesh/particle_and_bond_properties.txt','w');
fprintf(fileID,'%12.8f\n',particle_diameter);
fprintf(fileID,'%12.0f\n',num_particles);
fprintf(fileID,'%12.0f\n',num_bonds);

if FSI == 1
    fprintf(fileID,'%12.0f\n',num_boundary_particles);
    fprintf(fileID,'%12.0f\n',num_boundary_bonds);
end

fprintf(fileID,'%12.0f\n',num_clamped_particles);
fprintf(fileID,'%12.0f\n',resolution);
fclose(fileID);

%%%% File with particle positions
fileID = fopen('Mesh/particle_positions.txt','w');
fprintf(fileID,'%12.8f\n',particlePositionX);
fprintf(fileID,'%12.8f\n',particlePositionY);
fclose(fileID);

%%%% File with bond connections
fileID = fopen('Mesh/bond_IDs.txt','w');
fprintf(fileID,'%12.0f\n',b2p_ID(:,1));
fprintf(fileID,'%12.0f\n',b2p_ID(:,2));
fclose(fileID);

if FSI == 1
    %%%% File with boundary particle IDs (indexed with particle position file)
    fileID = fopen('Mesh/boundary_particle_IDs.txt','w');
    fprintf(fileID,'%12.0f\n',particle_boundary);
    fclose(fileID);

    %%%% File with boundary bond connections
    fileID = fopen('Mesh/boundary_bond_IDs.txt','w');
    fprintf(fileID,'%12.0f\n',b2p_ID_boundary_sequence(:,1));
    fprintf(fileID,'%12.0f\n',b2p_ID_boundary_sequence(:,2));
    fclose(fileID);
end

%%%% File with clamped particle IDs
fileID = fopen('Mesh/clamped_particle_IDs.txt','w');
fprintf(fileID,'%12.0f\n',clamped_particles);
fclose(fileID);

end

   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Post processing %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if output_figures == 1
    
    if beam == 1

        beam_figure_output( particlePositionX, particlePositionY, num_particles, particle_diameter, num_boundary_particles, particle_boundary, height, width, start_x, start_y, x_desired, y_desired, num_bonds, b2p_ID, b2p_ID_boundary_sequence, num_boundary_bonds, num_clamped_particles, clamped_particles) 

    end

    if valve == 1

        valve_figure_output( particlePositionX, particlePositionY, num_particles, base_num_particles, leaflet_num_particles, base_particle_diameter, leaflet_particle_diameter, num_boundary_particles, particle_boundary, start_x, start_y, leaflet_height, x_desired, y_desired, num_bonds, b2p_ID, b2p_ID_boundary_sequence, num_boundary_bonds, num_clamped_particles, clamped_particles ) %%num_boundary_particles, particle_boundary, height, width, start_x, start_y, x_desired, y_desired, num_bonds, b2p_ID)

    end

    %%%% Print figures to file
    print(1,'-dpdf','-r600','Mesh/Figures/Particle Discretisation.pdf')
    print(2,'-dpdf','-r600','Mesh/Figures/Particles and Bonds.pdf')
    print(3,'-dpdf','-r600','Mesh/Figures/Boundary Particles.pdf')
    print(4,'-dpdf','-r600','Mesh/Figures/Boundary Bonds.pdf')
    print(5,'-dpdf','-r600','Mesh/Figures/Boundary Conditions.pdf')

end


