function out = trapz4Dto1D(F,x,y,z)
%TRAPZ4DTO1D Integrate 4D grid over the first three axes
%   Integration of 4D (x,y,z,f) meshgrid along x,y,z
    out = trapz(z,trapz(y,trapz(x,F,1),2),3);
    out = squeeze(out);
end