function [q]=getq(N, F,qmin,L)
% Get array that describes the distributed load (assumes linear
% distribution
% Inputs:
%   L - length of the beam
%   N - number of elements
%   F - total  force of the 
%   qmin - the minimum value of the distributed load
% Outputs:
%   q- the distributed load array
%   ----------------------------
%   calculate force in the triangular part of the distribution
    f=F-qmin*L;
    %calculate the max value of the distributed load
    qmax=2*f/L+qmin;
    %create the array of values for the distributed load
    q=linspace(qmax,qmin,N+1);
end