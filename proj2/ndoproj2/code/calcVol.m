function [vol,g]= calcVol(R,L,N)
% calculate the volume and the volume gradient for the spar
% Inputs:
%   R- the design variable array
%   L- the length of the spar
%   N- the number of elements
% Outputs:
%   vol - the volume of the spar
%   g - the gradient fo the volume calculated using the complex step method
%   ----------------------------

    function [vol]=subs(r)
        %calculate the volume and the volume gradient for the spar
% Inputs:
%   r- the design variable array
% Outputs:
%   vol - the volume of the spar
%   ----------------------------
    [rin,rout]=getR(r);
    vol=0;
    for i=1:length(rin)
        vol=vol+(pi*rout(i)^2-pi*rin(i)^2)*L/(N);
    end
    end
 %set up complex step
    h=1e-30;
    g=zeros(size(R));
    %compute gradient 
    for j=1:length(g)
        comp=R;
    comp(j)=R(j)+complex(0,h);
    g(j)=imag(subs(comp))/h;
    end
%compute volume
vol=subs(R);
%g=imag(subs(comp))/h;
end