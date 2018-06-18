function [r1, r2, r3, r4] = c_umu_randtetra(L, E, P)
%
% [r1, r2, r3, r4] = c_umu_randtetra(L, E, P);
%
% Calculates random SC positions with characteristic parameters L, E, P.
% The origin of coordinates is at the centre of mass of the tetrahedron.
% The tetrahedron is rotated randomly in space.
%
% The theory is describe in the file Tetraeqs.tex
%
% L           Average distance between SC (km)
% E           Elongation
% P           Planarity
% r1..r4      Row vector (col 1-3 are SC pos in km)

     if E > 1  
        error(sprintf('E = %g, must be less than 1.', E));
     end 
     if P > 1  
        error(sprintf('P = %g, must be less than 1.', P));
     end 

     if E < 0  
        error(sprintf('E = %g, must be greater than 0.', E));
     end 
     if P < 0  
        error(sprintf('P = %g, must be greater than 0.', P));
     end 


     theta = acos(rand(1,1));  % Random angles to give uniform
     phi = 2*pi*rand(1,1);     % probability on the surface of
     psi = 2*pi*rand(1,1);     % the unit sphere

     st = sin(theta);
     ct = cos(theta);
     sf = sin(phi);
     cf = cos(phi);
     sp = sin(psi);
     cp = cos(psi);

     u = L/2;                % Square root of the eigenvalues
     v = u*(1-E);            % of the volumetric tensor
     w = v*(1-P);

     alpha = u*[st*cf;          st*sf;          ct;     0];
     beta  = v*[ct*cf*cp-sf*sp; ct*sf*cp+cf*sp; -st*cp; 0];
     gamma = w*[ct*cf*sp+sf*cp; ct*sf*sp-cf*cp; -st*sp; 0];

     r = [1 1 1 1; -1 -1 1 1; 1 -1 -1 1; -1 1 -1 1]*[alpha beta gamma];

     theta = acos(rand(1,1));        % Random rotation
     phi = 2*pi*rand(1,1);
     psi = 2*pi*rand(1,1);
     [R] = c_umu_rotmatr(theta, phi, psi);
     r = R*r';
     r1 = r(:, 1)';
     r2 = r(:, 2)';
     r3 = r(:, 3)';
     r4 = r(:, 4)';

