function [c,ceq,j,jeq]=getC(R,L,m,N,E,sigmaUlt)
% calculate the nonlinear inequality constraint and the nonlinear
%   inequality constraint gradient for the spar
% Inputs:
%   R- the design variable array
%   L- the length of the spar
%   N- the number of elements
%   m-the mass of the plane
%   E- he young's modulus of the spar material
%   sigmaUlt - the ultimate tensile strength of the material
% Outputs:
%   c - the normalized nonlinear inequality constraint
%   j - the jacobian for the nonlinear inequality constraint calculated
%       using the complex step method
%   ----------------------------
    function [frac]=subs(r)
        % calculate the nonlinear inequality constraint 
% Inputs:
%   r- the design variable array
% Outputs:
%   c - the normalized nonlinear inequality constraint
%   ----------------------------
    %get the distributed load array
    q=getq(N, 2.5*m*9.81,0,L);
    %get rin, rout
    [rin,rout]=getR(r);
    %get the second moment of Inertia array
    Iyy=getI(r);
    %get the beam displacements
    [u] = CalcBeamDisplacement(L, E, Iyy, q, N);
    %get element stresses
    [sigma] = CalcBeamStress(L, E, rout, u, N);
    %create c 
    frac=sigma/sigmaUlt;

    end
    %set up complex step
    h=1e-30;
    g=zeros([2*(N+1),N+1]);
    %calculate jacobian
    for i=1:2*(N+1)
    comp=R;
    comp(i)=R(i)+complex(0,h);
    g(i,:)= imag(subs(comp))/h;
    end
    frac=subs(R);
    %calculate outputs
    c=frac-ones(size(frac));
    ceq=[];
    jeq=[];
    j=g;


end
