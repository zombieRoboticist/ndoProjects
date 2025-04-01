function [rin,rout]=getR(R)
%Splits a design variable into rin and rout
% Inputs:
%   R- the design variable array
% Outputs:
%   rin - the array of inner radii
%   rout - the array of outer radii
%   ----------------------------
%
l=max(size(R))/2;
    rout=zeros([l,1]);
    rin=zeros([l,1]);
    for i=1:l
        rin(i)=R(2*i-1);
        rout(i)=R(2*i);
    end
end