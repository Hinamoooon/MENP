function [Cp,CT,Cm,CQe,CQm,Csum] = toroidalME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z)
%toroidME Multipole Expansion with toroidal dipole
% Multipole expansion is calculated from current density distribution in
% Cartesian coordinates. Toroidal dipole is introduced.
% Note that higher order terms, including mean-square radii multipoles
% (m(1),T(1),...), are neglected in scattering cross sections.
% [ref.1,Table1]
% For Qe and Qm, expressions in [ref.2,Table1] is adopted, that is,
% high-order toroidal (quadrupole) moment is included in Qe.
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
% CT: toroidal dipole
% Cm: magnetic dipole
% CQe: electric quadrupole
% CQm: magnetic quadrupole
%
% References
% 1. "Optical Anapoles: Concepts and Applications"
%    (http://doi.wiley.com/10.1002/adom.201801350)
% 2. "An electromagnetic multipole expansion beyond the long-wavelength approximation"
%    (http://dx.doi.org/10.1016/j.optcom.2017.08.064)
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
    % constant for scattering cross section
    const = k.^4/(6*pi*eps0^2*1);  % E0 = 1
    
    % scalar product
    rJ = x4d.*Jx + y4d.*Jy + z4d.*Jz;  % product r,J(r)
    rr = x4d.*x4d + y4d.*y4d + z4d.*z4d;  % product r,r
    
    % cross product r x J = (ry*Jz-rz*Jy, rz*Jx-rx*Jz, rx*Jy-ry*Jx)
    rxJx = (y4d.*Jz - z4d.*Jy);
    rxJy = (z4d.*Jx - x4d.*Jz);
    rxJz = (x4d.*Jy - y4d.*Jx);
    
    %% calculate multipole moments and cross sections
    % calculate electric dipole p,
    px = -1./(1i*omega).*(trapz4Dto1D(Jx,x,y,z));
    py = -1./(1i*omega).*(trapz4Dto1D(Jy,x,y,z));
    pz = -1./(1i*omega).*(trapz4Dto1D(Jz,x,y,z));
    norm2_p = px.*conj(px)+py.*conj(py)+pz.*conj(pz);
    Cp = const.*norm2_p;
    
    % calculate toroidal dipole T,
    dTx = rJ.*x4d-2*rr.*Jx;
    dTy = rJ.*y4d-2*rr.*Jy;
    dTz = rJ.*z4d-2*rr.*Jz;
    Tx = 1/(10*c)*trapz4Dto1D(dTx,x,y,z);
    Ty = 1/(10*c)*trapz4Dto1D(dTy,x,y,z);
    Tz = 1/(10*c)*trapz4Dto1D(dTz,x,y,z);
    norm2_T = Tx.*conj(Tx)+Ty.*conj(Ty)+Tz.*conj(Tz);
    CT = const.*abs((1i*k).^2.*norm2_T);
    
    % calculate net electric dipole moment (p+ikT) pT,
    pTx = px+1i*k.*Tx;
    pTy = py+1i*k.*Ty;
    pTz = pz+1i*k.*Tz;
    norm2_pT = pTx.*conj(pTx)+pTy.*conj(pTy)+pTz.*conj(pTz);
    CpT = const.*norm2_pT;
    
    % calculate magnetic dipole m,
    % dm:=r x J(r)
    dmx = rxJx;
    dmy = rxJy;
    dmz = rxJz;
    mx = 1/2*trapz4Dto1D(dmx,x,y,z);
    my = 1/2*trapz4Dto1D(dmy,x,y,z);
    mz = 1/2*trapz4Dto1D(dmz,x,y,z);
    norm2_m = mx.*conj(mx)+my.*conj(my)+mz.*conj(mz);
    Cm = const.*norm2_m/c^2;

    % calculate electric quadrupole Qe
    dQe1xx = 3*2*x4d.*Jx - 2*rJ;
    dQe1xy = 3*(y4d.*Jx+x4d.*Jy);
    dQe1xz = 3*(z4d.*Jx+x4d.*Jz);
    dQe1yy = 3*2*y4d.*Jy - 2.*rJ;
    dQe1yx = 3*(x4d.*Jy+y4d.*Jx);
    dQe1yz = 3*(z4d.*Jy+y4d.*Jz);
    dQe1zz = 3*2*z4d.*Jz - 2*rJ;
    dQe1zx = 3*(x4d.*Jz+z4d.*Jx);
    dQe1zy = 3*(y4d.*Jz+z4d.*Jy);
    dQe2xx = 4*x4d.*x4d.*rJ - 5*rr*2.*x4d.*Jx + 2*rr.*rJ;
    dQe2xy = 4*x4d.*y4d.*rJ - 5*rr.*(x4d.*Jy+y4d.*Jx);
    dQe2xz = 4*x4d.*z4d.*rJ - 5*rr.*(x4d.*Jz+z4d.*Jx);
    dQe2yy = 4*y4d.*y4d.*rJ - 5*rr*2.*y4d.*Jy + 2*rr.*rJ;
    dQe2yx = 4*y4d.*x4d.*rJ - 5*rr.*(y4d.*Jx+x4d.*Jy);
    dQe2yz = 4*y4d.*z4d.*rJ - 5*rr.*(y4d.*Jz+z4d.*Jy);
    dQe2zz = 4*z4d.*z4d.*rJ - 5*rr*2.*z4d.*Jz + 2*rr.*rJ;
    dQe2zx = 4*z4d.*x4d.*rJ - 5*rr.*(z4d.*Jx+x4d.*Jz);
    dQe2zy = 4*z4d.*y4d.*rJ - 5*rr.*(z4d.*Jy+y4d.*Jz);
    Qexx = -1./(1i*omega).*(trapz4Dto1D(dQe1xx,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2xx,x,y,z));
    Qexy = -1./(1i*omega).*(trapz4Dto1D(dQe1xy,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2xy,x,y,z));
    Qexz = -1./(1i*omega).*(trapz4Dto1D(dQe1xz,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2xz,x,y,z));
    Qeyy = -1./(1i*omega).*(trapz4Dto1D(dQe1yy,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2yy,x,y,z));
    Qeyx = -1./(1i*omega).*(trapz4Dto1D(dQe1yx,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2yx,x,y,z));
    Qeyz = -1./(1i*omega).*(trapz4Dto1D(dQe1yz,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2yz,x,y,z));
    Qezz = -1./(1i*omega).*(trapz4Dto1D(dQe1zz,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2zz,x,y,z));
    Qezx = -1./(1i*omega).*(trapz4Dto1D(dQe1zx,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2zx,x,y,z));
    Qezy = -1./(1i*omega).*(trapz4Dto1D(dQe1zy,x,y,z)+ ...
           k.^2/14.*trapz4Dto1D(dQe2zy,x,y,z));
    norm2_Qe = Qexx.*conj(Qexx)+Qexy.*conj(Qexy)+Qexz.*conj(Qexz)+ ...
               Qeyy.*conj(Qeyy)+Qeyx.*conj(Qeyx)+Qeyz.*conj(Qeyz)+ ...
               Qezz.*conj(Qezz)+Qezx.*conj(Qezx)+Qezy.*conj(Qezy);
    CQe = const/120.*k.^2.*norm2_Qe;

    % calculate magnetic quadrupole Qm
    dQmxx = 2*x4d.*rxJx;
    dQmxy = x4d.*rxJy+y4d.*rxJx;
    dQmxz = x4d.*rxJz+x4d.*rxJz;
    dQmyy = 2*y4d.*rxJy;
    dQmyx = y4d.*rxJx+x4d.*rxJy;
    dQmyz = y4d.*rxJz+z4d.*rxJy;
    dQmzz = 2*z4d.*rxJz;
    dQmzx = z4d.*rxJx+x4d.*rxJz;
    dQmzy = z4d.*rxJy+y4d.*rxJz;
    Qmxx = trapz4Dto1D(dQmxx,x,y,z);
    Qmxy = trapz4Dto1D(dQmxy,x,y,z);
    Qmxz = trapz4Dto1D(dQmxz,x,y,z);
    Qmyy = trapz4Dto1D(dQmyy,x,y,z);
    Qmyx = trapz4Dto1D(dQmyx,x,y,z);
    Qmyz = trapz4Dto1D(dQmyz,x,y,z);
    Qmzz = trapz4Dto1D(dQmzz,x,y,z);
    Qmzx = trapz4Dto1D(dQmzx,x,y,z);
    Qmzy = trapz4Dto1D(dQmzy,x,y,z);
    norm2_Qm = Qmxx.*conj(Qmxx)+Qmxy.*conj(Qmxy)+Qmxz.*conj(Qmxz)+ ...
               Qmyy.*conj(Qmyy)+Qmyx.*conj(Qmyx)+Qmyz.*conj(Qmyz)+ ...
               Qmzz.*conj(Qmzz)+Qmzx.*conj(Qmzx)+Qmzy.*conj(Qmzy);
    CQm = const./120.*(k/c).^2.*norm2_Qm;
    
    % sum of all multipoles
    Csum = CpT + Cm + CQe + CQm;

end