function [phidot]=dynamics(tau,phi0,r2, alpha1,omega)
%calculate the derivative of the dynamics
% Inputs:
%   tau - the nondimensional time to evaluate the dynamics at
%   phi0 - the angular velocity and position of the car
%   r2 - the distance between the cars center of gravity and piviot point
%   alpha1 - the change in angle with respect to tau
%   omega - the angular velocity of the ride
% Outputs:
%   phidot - the angular velocity and acceleration of the car
%   ----------------------------

%set up constants
q0=20;
r1=4.3; 
alpha0=.036;
%calculate coefficients
alpha= alpha0-alpha1*cos(tau);
beta=3*alpha1*sin(tau);
gamma=sqrt(9.81/r2)/(3*omega);
e=r1/(9*r2);

%calculate dynamics
phidot=[phi0(2);-1*((e-gamma*gamma*alpha)*sin(phi0(1))+gamma*gamma*beta*cos(phi0(1))+gamma*phi0(2)/q0)];