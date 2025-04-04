function [sigma] = obj(X, numRot)
%calculate the standard deviation of the dynamics
% Inputs:
%   X0 - (r2) the distance between the cars center of gravity and piviot point
%   X1 - (alpha1) the change in angle with respect to tau
%   X2 - (omega) the angular velocity of the ride
%   numRot - the number of rotations to evaluate the ride at
% Outputs:
%   sigma - the standard deviation of the dynamics
%   ----------------------------


%set ode 45 vars
%tspan from 0 to 2pi*numRot because it goes in numRot full circles 
tspan=[0,numRot*2*pi];
%x0 -pi/2 because the ride at rest should point down hill, 0 bc at rest
x0=[-pi/2,0];
odefunc=@(t,y) dynamics(t,y,X(1),X(2),X(3));
%run ode45

[t,y]=ode45(odefunc,tspan,x0);
%plot(t,y(:,2))
%hold on
%plot([0,tspan(2)],[mean(y(:,2)),mean(y(:,2))])
%plot(t,y(:,1))
%hold off
%ouput standard deviation
sigma=3*X(3)*sqrt(trapz(t,((y(:,2)-trapz(t,(y(:,2)))/tspan(2)).^2))/(tspan(2)));
%sigma=3*X(3)*sqrt(trapz(t,((y(:,2)-trapz(t,(y(:,2)))/max(size(y))).^2))/(tspan(2)));

%

