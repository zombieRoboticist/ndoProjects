function [Iyy]=getI(R)
    %Calculates Iyy for each element in the model
% Inputs:
%   R- the design variable array
% Outputs:
%   Iyy - the array of the second moment of inertia
%   ----------------------------
%get rin, rout from design variable array
    [rin,rout]=getR(R);
    Iyy=zeros(length(rin));
%calculate Iyy for each element
    for i=1:length(rin)
        Iyy(i)=pi*(rout(i)^4-rin(i)^4)/4;
    end
end