function [Cp,Cm,CQe,CQm,Csum] = exactME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z)
%EXACTME Exact Multipole Expansion
% Exact multipole expansion (beyond long-wavelength approximation)
% is calculated from current density distribution.[ref.1,Table2]
%
% Input properties
% x,y,z: 1D arrays which defines position vector r=(x,y,z) toward each
%        meshgrid point to be calculated. Units are [m].
% Ex,Ey,Ez: Complex electric fields E(r,f)=E(x,y,z,f)=(Ex,Ey,Ez).
%           Exported data from EM simulation software can be used.
% n_x,n_y,n_z: Refractive indices at the corresponding meshgrid.
% 
% Output properties
% Cp,Cm,CQe,CQm: Scattering cross section from each of multipole moment.
% Cp: electric dipole
% Cm: magnetic dipole
% CQe: electric quadrupole
% CQm: magnetic quadrupole
%
% References
% 1. "An electromagnetic multipole expansion beyond the long-wavelength approximation"
%    (http://dx.doi.org/10.1016/j.optcom.2017.08.064)
%
% MENP (Multipole Expansion for NanoPhotonics)
% T. Hinamoto (Kobe University, Japan)

    %% prepare
    PhysConst;

    [x4d,y4d,z4d,f4d] = ndgrid(x,y,z,f);
    
    k4d = 2*pi*f4d/c;
    omega = 2*pi*f;
    k = omega/c;

    %% get current density
    [Jx,Jy,Jz] = E2J(Ex,Ey,Ez,n_x,n_y,n_z,f4d);

    %% calculate often used values
    % constant for scattering cross section
    const = k.^4/(6*pi*eps0^2*1);  % E0 = 1
    
    % scalar product
    rJ = x4d.*Jx + y4d.*Jy + z4d.*Jz;  % product r,J(r)
    rr = x4d.*x4d + y4d.*y4d + z4d.*z4d;  % product r,r
    r = sqrt(rr);  % norm(r)
    
    % cross product r x J = (ry*Jz-rz*Jy, rz*Jx-rx*Jz, rx*Jy-ry*Jx)
    rxJx = (y4d.*Jz - z4d.*Jy);
    rxJy = (z4d.*Jx - x4d.*Jz);
    rxJz = (x4d.*Jy - y4d.*Jx);

    % spherical bessel functions
    sbj0 = sqrt(pi./(2*k4d.*r)).*besselj(0+1/2,k4d.*r);  %j0
    sbj1 = sqrt(pi./(2*k4d.*r)).*besselj(1+1/2,k4d.*r);  %j1
    sbj2 = sqrt(pi./(2*k4d.*r)).*besselj(2+1/2,k4d.*r);  %j2
    sbj3 = sqrt(pi./(2*k4d.*r)).*besselj(3+1/2,k4d.*r);  %j3
    
    %% calculate multipole moments and cross sections
    % calculate electric dipole p,
    dpx = (3*rJ.*x4d-rr.*Jx).*sbj2./(k4d.*r).^2;
    dpy = (3*rJ.*y4d-rr.*Jy).*sbj2./(k4d.*r).^2;
    dpz = (3*rJ.*z4d-rr.*Jz).*sbj2./(k4d.*r).^2;
    px = -1./(1i*omega).*(trapz4Dto1D(Jx.*sbj0,x,y,z)+ ...
         k.^2/2.*trapz4Dto1D(dpx,x,y,z));
    py = -1./(1i*omega).*(trapz4Dto1D(Jy.*sbj0,x,y,z)+ ...
         k.^2/2.*trapz4Dto1D(dpy,x,y,z));
    pz = -1./(1i*omega).*(trapz4Dto1D(Jz.*sbj0,x,y,z)+ ...
         k.^2/2.*trapz4Dto1D(dpz,x,y,z));
    norm2_p = px.*conj(px)+py.*conj(py)+pz.*conj(pz);
    Cp = const.*norm2_p;

    % calculate magnetic dipole m,
    dmx = rxJx.*sbj1./(k4d.*r);
    dmy = rxJy.*sbj1./(k4d.*r);
    dmz = rxJz.*sbj1./(k4d.*r);
    mx = 3/2*trapz4Dto1D(dmx,x,y,z);
    my = 3/2*trapz4Dto1D(dmy,x,y,z);
    mz = 3/2*trapz4Dto1D(dmz,x,y,z);
    norm2_m = mx.*conj(mx)+my.*conj(my)+mz.*conj(mz);
    Cm = const.*norm2_m/c^2;

    % calculate electric quadrupole Qe
    dQe1xx = (3*2*x4d.*Jx - 2*rJ).*sbj1./(k4d.*r);
    dQe1xy = (3*(y4d.*Jx+x4d.*Jy)).*sbj1./(k4d.*r);
    dQe1xz = (3*(z4d.*Jx+x4d.*Jz)).*sbj1./(k4d.*r);
    dQe1yy = (3*2*y4d.*Jy - 2.*rJ).*sbj1./(k4d.*r);
    dQe1yx = (3*(x4d.*Jy+y4d.*Jx)).*sbj1./(k4d.*r);
    dQe1yz = (3*(z4d.*Jy+y4d.*Jz)).*sbj1./(k4d.*r);
    dQe1zz = (3*2*z4d.*Jz - 2*rJ).*sbj1./(k4d.*r);
    dQe1zx = (3*(x4d.*Jz+z4d.*Jx)).*sbj1./(k4d.*r);
    dQe1zy = (3*(y4d.*Jz+z4d.*Jy)).*sbj1./(k4d.*r);
    dQe2xx = (5*x4d.*x4d.*rJ - rr*2.*x4d.*Jx - rr.*rJ).*sbj3./(k4d.*r).^3;
    dQe2xy = (5*x4d.*y4d.*rJ - rr.*(x4d.*Jy+y4d.*Jx)).*sbj3./(k4d.*r).^3;
    dQe2xz = (5*x4d.*z4d.*rJ - rr.*(x4d.*Jz+z4d.*Jx)).*sbj3./(k4d.*r).^3;
    dQe2yy = (5*y4d.*y4d.*rJ - rr*2.*y4d.*Jy - rr.*rJ).*sbj3./(k4d.*r).^3;
    dQe2yx = (5*y4d.*x4d.*rJ - rr.*(y4d.*Jx+x4d.*Jy)).*sbj3./(k4d.*r).^3;
    dQe2yz = (5*y4d.*z4d.*rJ - rr.*(y4d.*Jz+z4d.*Jy)).*sbj3./(k4d.*r).^3;
    dQe2zz = (5*z4d.*z4d.*rJ - rr*2.*z4d.*Jz - rr.*rJ).*sbj3./(k4d.*r).^3;
    dQe2zx = (5*z4d.*x4d.*rJ - rr.*(z4d.*Jx+x4d.*Jz)).*sbj3./(k4d.*r).^3;
    dQe2zy = (5*z4d.*y4d.*rJ - rr.*(z4d.*Jy+y4d.*Jz)).*sbj3./(k4d.*r).^3;
    Qexx = -3./(1i*omega).*(trapz4Dto1D(dQe1xx,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2xx,x,y,z));
    Qexy = -3./(1i*omega).*(trapz4Dto1D(dQe1xy,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2xy,x,y,z));
    Qexz = -3./(1i*omega).*(trapz4Dto1D(dQe1xz,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2xz,x,y,z));
    Qeyy = -3./(1i*omega).*(trapz4Dto1D(dQe1yy,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2yy,x,y,z));
    Qeyx = -3./(1i*omega).*(trapz4Dto1D(dQe1yx,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2yx,x,y,z));
    Qeyz = -3./(1i*omega).*(trapz4Dto1D(dQe1yz,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2yz,x,y,z));
    Qezz = -3./(1i*omega).*(trapz4Dto1D(dQe1zz,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2zz,x,y,z));
    Qezx = -3./(1i*omega).*(trapz4Dto1D(dQe1zx,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2zx,x,y,z));
    Qezy = -3./(1i*omega).*(trapz4Dto1D(dQe1zy,x,y,z)+ ...
           2*k.^2.*trapz4Dto1D(dQe2zy,x,y,z));
    norm2_Qe = Qexx.*conj(Qexx)+Qexy.*conj(Qexy)+Qexz.*conj(Qexz)+ ...
               Qeyy.*conj(Qeyy)+Qeyx.*conj(Qeyx)+Qeyz.*conj(Qeyz)+ ...
               Qezz.*conj(Qezz)+Qezx.*conj(Qezx)+Qezy.*conj(Qezy);
    CQe = const/120.*k.^2.*norm2_Qe;

    % calculate magnetic quadrupole Qm
    dQmxx = (2*x4d.*rxJx).*sbj2./(k4d.*r).^2;
    dQmxy = (x4d.*rxJy+y4d.*rxJx).*sbj2./(k4d.*r).^2;
    dQmxz = (x4d.*rxJz+x4d.*rxJz).*sbj2./(k4d.*r).^2;
    dQmyy = (2*y4d.*rxJy).*sbj2./(k4d.*r).^2;
    dQmyx = (y4d.*rxJx+x4d.*rxJy).*sbj2./(k4d.*r).^2;
    dQmyz = (y4d.*rxJz+z4d.*rxJy).*sbj2./(k4d.*r).^2;
    dQmzz = (2*z4d.*rxJz).*sbj2./(k4d.*r).^2;
    dQmzx = (z4d.*rxJx+x4d.*rxJz).*sbj2./(k4d.*r).^2;
    dQmzy = (z4d.*rxJy+y4d.*rxJz).*sbj2./(k4d.*r).^2;
    Qmxx = 15*trapz4Dto1D(dQmxx,x,y,z);
    Qmxy = 15*trapz4Dto1D(dQmxy,x,y,z);
    Qmxz = 15*trapz4Dto1D(dQmxz,x,y,z);
    Qmyy = 15*trapz4Dto1D(dQmyy,x,y,z);
    Qmyx = 15*trapz4Dto1D(dQmyx,x,y,z);
    Qmyz = 15*trapz4Dto1D(dQmyz,x,y,z);
    Qmzz = 15*trapz4Dto1D(dQmzz,x,y,z);
    Qmzx = 15*trapz4Dto1D(dQmzx,x,y,z);
    Qmzy = 15*trapz4Dto1D(dQmzy,x,y,z);
    norm2_Qm = Qmxx.*conj(Qmxx)+Qmxy.*conj(Qmxy)+Qmxz.*conj(Qmxz)+ ...
               Qmyy.*conj(Qmyy)+Qmyx.*conj(Qmyx)+Qmyz.*conj(Qmyz)+ ...
               Qmzz.*conj(Qmzz)+Qmzx.*conj(Qmzx)+Qmzy.*conj(Qmzy);
    CQm = const./120.*(k/c).^2.*norm2_Qm;
    
    % sum of all multipoles
    Csum = Cp + Cm + CQe + CQm;
    
end