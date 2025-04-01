function [R]=putR(rin, rout)
%Integrate rin and rout into a R array
% Inputs:
%   rin - the array of inner radii
%   rout - the array of outer radii
% Outputs:
%   R- the design variable array
%   ----------------------------
    R=zeros([length(rin)*2,1]);
    for i=1:length(rin)
        R(2*i-1)=rin(i);
        R(2*i)=rout(i);
    end
end
