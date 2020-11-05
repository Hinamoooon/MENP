function [Jx,Jy,Jz] = E2J(Ex,Ey,Ez,n_x,n_y,n_z,f4d)
%E2J calculates current density J=(Jx,Jy,Jz) from electric fields
%   Input
%   f: 4D meshgrid of frequency in [Hz].
%   Ex,Ey,Ez: Complex electric fields at the each grid points (4D matrix).
%   n_x,n_y,n_z: Refractive indices at the grid points (4D matrix).
%
% MENP (Multipole Expansion for NanoPhotonics)
% T. Hinamoto (Kobe University, Japan)
    
    PhysConst;
    
    % J(r)=-i*omega*eps0*(n^2-1)E(r), where omega = 2*pi*f
    Jx = - 1i * 2*pi*f4d * eps0 .* (n_x.^2 - 1) .* Ex;
    Jy = - 1i * 2*pi*f4d * eps0 .* (n_y.^2 - 1) .* Ey;
    Jz = - 1i * 2*pi*f4d * eps0 .* (n_z.^2 - 1) .* Ez;

end