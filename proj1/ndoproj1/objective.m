function [f]=objective(a,Nx, Ny, L, kappa, Ttop,Tbot)
% Solves for the heat flux using the CalcFlux Function, and returns the 
% inverse heat flux per unit length from the water to the air
% Inputs:
%   L - length of domain in x direction
%  a - the coefficient vector for the shape of the heat exchanger
%   Nx - number of elements along the x direction
%   Ny - number of elements along the y direction
%   kappa - thermal conductivity
%   Ttop - the ambient air temperature along the top of the domain
%   Tbot - the fluid temperature along the bottom of the domain
% Outputs:
%   f -the inverse of heat flux per unit length out of top/bottom boundary
% Notes:
%   The domain has straight sides on the left and right and bottom.  The
%   boundary conditions are periodic on the left and right.  Prescribed
%   (i.e. Dirichlet) boundary conditions are applied along the top and
%   bottom.  The flux is computed along the bottom boundary.
%--------------------------------------------------------------------------
%create the position array 
X = [0:L/Nx:L].';
%create the height array that defines the upper part of the heat exchanger
h=calcHeight(a,L,X);
%calculate the heat flux/unit length using the CalcFlux function
[flux,T,dTdx,xy] = CalcFlux(L, h, Nx, Ny, kappa, Ttop, Tbot);
%calculate the inverse of the flux and set it to f
f=1/flux;
end