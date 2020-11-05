function [arg_px,arg_py,arg_pz,arg_ikTx,arg_ikTy,arg_ikTz] = toroidalME_phase(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z)
%toroidalME_phase Phase of electric and toroidal dipole moments
% Multipole expansion is calculated from current density distribution in
% Cartesian coordinates. Toroidal dipole is introduced.
% [ref.1,Table1]
%
% Input properties
% x,y,z: 1D arrays which defines position vector r=(x,y,z) toward each
%        meshgrid point to be calculated. Units are [m].
% Ex,Ey,Ez: Complex electric fields E(r,f)=E(x,y,z,f)=(Ex,Ey,Ez).
%           Exported data from EM simulation software can be used.
% n_x,n_y,n_z: Refractive indices at the corresponding meshgrid.
% 
% Output properties
% arg_p, arg_ikT: phases of  electric and toroidal (as -ikT) dipole moments
% arg_px, arg_py, arg_pz: electric dipole
% arg_ikTx, arg_ikTy, arg_ikTz: toroidal dipole
%
% References
% 1. "Optical Anapoles: Concepts and Applications"
%    (http://doi.wiley.com/10.1002/adom.201801350)
%
% MENP (Multipole Expansion for NanoPhotonics)
% T. Hinamoto (Kobe University, Japan)

    %% prepare
    PhysConst;

    [x4d,y4d,z4d,f4d] = ndgrid(x,y,z,f);
    
    omega = 2*pi*f;
    k = omega/c;

    %% get current density
    [Jx,Jy,Jz] = E2J(Ex,Ey,Ez,n_x,n_y,n_z,f4d);

    %% calculate often used values
    % scalar product
    rJ = x4d.*Jx + y4d.*Jy + z4d.*Jz;  % product r,J(r)
    rr = x4d.*x4d + y4d.*y4d + z4d.*z4d;  % product r,r
    
    %% calculate multipole moments and cross sections
    % calculate electric dipole p,
    px = -1./(1i*omega).*(trapz4Dto1D(Jx,x,y,z));
    py = -1./(1i*omega).*(trapz4Dto1D(Jy,x,y,z));
    pz = -1./(1i*omega).*(trapz4Dto1D(Jz,x,y,z));
    arg_px = angle(px);
    arg_py = angle(py);
    arg_pz = angle(pz);
%     arg_px = atan(real(px)./imag(px));
%     arg_py = atan(real(py)./imag(py));
%     arg_pz = atan(real(pz)./imag(pz));
    
    % calculate toroidal dipole T,
    dTx = rJ.*x4d-2*rr.*Jx;
    dTy = rJ.*y4d-2*rr.*Jy;
    dTz = rJ.*z4d-2*rr.*Jz;
    Tx = 1/(10*c)*trapz4Dto1D(dTx,x,y,z);
    Ty = 1/(10*c)*trapz4Dto1D(dTy,x,y,z);
    Tz = 1/(10*c)*trapz4Dto1D(dTz,x,y,z);
    ikTx = -1i*k.*Tx;
    ikTy = -1i*k.*Ty;
    ikTz = -1i*k.*Tz;
    arg_ikTx = angle(ikTx);
    arg_ikTy = angle(ikTy);
    arg_ikTz = angle(ikTz);
%     arg_ikTx = atan(real(ikTx)./imag(ikTx));
%     arg_ikTy = atan(real(ikTy)./imag(ikTy));
%     arg_ikTz = atan(real(ikTz)./imag(ikTz));

end