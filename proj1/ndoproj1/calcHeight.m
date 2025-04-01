function calcHeight=calcHeight(a,L,X)
% Creates the height of the 
% Inputs:
%   L - length of domain in x direction
%  a - the coefficient vector for the shape of the heat exchanger
%  X - the position array
% Outputs:
%   h - the height array that defines the shape of the upper part of the heat exchanger
%--------------------------------------------------------------------------
%initialize the height array and set it to a0
h=ones(size(X))*a(1);
%loops through every position in X to calculate the height of the 
for i = 2:max(size(X))
        %loops through the coefficients in a to calculate the h at a 
        % specific position X(i)
        for j=1:max(size(a))
            %at each position X(i) and coefficient a(j) add the value of
            %the sinusoidal function a*sin((2*pi*(j-1)*X)/L)
            h(i) = h(i)+a(j)*sin(2*pi*(j-1)*X(i)/L);
        end
end
%h(1)=a(1);
%return the value of the height array 
calcHeight=h;
end