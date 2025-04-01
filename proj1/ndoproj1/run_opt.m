function [a]=run_opt(Nx,Ny,a0)
% Sets up and runs the optimzation algorithm fmincon
% Inputs:
%   Nx - number of elements along the x direction
%   Ny - number of elements along the y direction
%   a0 - the initial guess for the coefficients
% Outputs:
%   a - the resulting optimized coefficients from the optimizer
%
% Notes:
%   The domain has straight sides on the left and right and bottom.  The
%   boundary conditions are periodic on the left and right.  Prescribed
%   (i.e. Dirichlet) boundary conditions are applied along the top and
%   bottom.  The flux is computed along the bottom boundary.
%--------------------------------------------------------------------------
%temperature of the water
Twater=90;
%temperature of the air
Tair=20;
%length of the repeating segment
L=.05;
%thermal conductivity of the heat exchanger
kappa=20;
%create the position array 
x = [0:L/Nx:L].';
%set fmincon options
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');




%create anonomus function
func = @(a) objective(a,Nx, Ny, L, kappa, Tair,Twater);

%define size constraints
nvar = size(a0,1); % number of design variables
Aineq = zeros(2*(Nx+1),nvar); 
bineq = zeros(2*(Nx+1),1);
for i = 1:(Nx+1)
  % first, the upper bound
  Aineq(i,1) = 1.0; % this coefficient corresponds to variable a_1
  for k = 2:nvar
    Aineq(i,k) = sin(2*pi*(k-1)*x(i)/L); % this coefficient corresponds to variable a_k
  end
  bineq(i) = 0.05; % the upper bound value
  % next, the lower bound; we use ptr to shift the index in Aineq and bineq
  ptr = Nx+1;
  Aineq(ptr+i,1) = -1.0; % note the negative here!!! fmincon expects inequality in a form A*x < b
  for k = 2:nvar
    Aineq(ptr+i,k) = -sin(2*pi*(k-1)*x(i)/L); % again, a negative
  end
  bineq(ptr+i) = -0.01; % negative lower bound
end

%call fmincon and store the output in a
a = fmincon(func,a0,Aineq,bineq,[],[],[],[],[],options);
end
