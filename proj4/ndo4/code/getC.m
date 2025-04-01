function [c,ceq,j,jeq]=getC(R,L,m,N,E,sigmaUlt,ngq)
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
    function [bob]=stressX(q,x,N,L,r,E)
    % calculate the nonlinear inequality constraint 
% Inputs:
%   r- the design variable array
%   L- the length of the spar
%   N- the number of elements
%   q-the nominal loding condition
%   E- he young's modulus of the spar material
%   x- the position along the spar
% Outputs:
%   bob - the uncertain spar loading
%   ----------------------------
    for i=1:N+1
        q(i)=q(i)+x(1)*cos(pi*(i-1)/2/(N))+x(2)*cos(3*pi*(i-1)/2/(N))+x(3)*cos(5*pi*(i-1)/2/(N))+x(4)*cos(7*pi*(i-1)/2/(N));
    end
    %get rin, rout
    [rin,rout]=getR(r);
    %get the second moment of Inertia array
    Iyy=getI(r);   
    %get the beam displacements
    [u] = CalcBeamDisplacement(L, E, Iyy, q, N);
    %get element stresses
    [bob] = CalcBeamStress(L, E, rout, u, N);
    bob;
end
    %get the distributed load array
    q=getq(N, 2.5*m*9.81,0,L);
    %gauss quadrature points
    sigf=[q(1)/10,q(1)/20,q(1)/30,q(1)/40];
    if ngq == 1
        xi = [0.0];
        wts = [1.77245]./sqrt(pi);% adjusted weights !
    elseif ngq ==2
        xi = [-0.707107; 0.707107];
        wts = [0.886227; 0.886227]./sqrt(pi);
    elseif ngq == 3
        xi = [ -1.22474487139; 0.0; 1.22474487139];
        wts = [0.295408975151; 1.1816359006; 0.295408975151]./ sqrt(pi); 
    elseif ngq == 4
        xi = [-1.65068; -0.524648; 0.524648; 1.65068];
        wts = [0.0813128; 0.804914; 0.804914; 0.0813128]./sqrt(pi);
    elseif ngq == 5
        xi = [-2.02018; -0.958572; 0.0; 0.958572; 2.02018];
        wts = [0.0199532; 0.393619; 0.945309; 0.393619; 0.0199532]./sqrt(pi);
    elseif ngq == 6
        xi = [-2.35061; -1.33585; -0.436077; 0.436077; 1.33585; 2.35061];
        wts = [0.00453001; 0.157067; 0.72463; 0.72463; 0.157067; 0.00453001]./sqrt(pi);
    end
    mean=0;
    mean2=0;
    %calculate mean, mean squared using gauss quadrature
    for i1=1:ngq
        p1=sqrt(2)*sigf(1)*xi(i1);
        for i2=1:ngq
            p2=sqrt(2)*sigf(2)*xi(i2);
            for i3=1:ngq
                p3=sqrt(2)*sigf(3)*xi(i3);
                for i4=1:ngq
                    p4=sqrt(2)*sigf(4)*xi(i4);
                    temp=stressX(q,[p1,p2,p3,p4],N,L,r,E);
                    mean=mean+wts(i1)*wts(i2)*wts(i3)*wts(i4)*temp;
                    mean2=mean2+wts(i1)*wts(i2)*wts(i3)*wts(i4)*temp.*temp;
                end
            end
        end
    end
    %calculate standard deviation
    sigma=mean+6*sqrt((mean2-mean.*mean));
    %create c 
    %plot mean +/- 6stddev
   % x=0:L/N:L;
    %plot(x,mean)
    %hold on
    %plot(x,mean+6*sqrt((mean2-mean.*mean)))
    %plot(x,mean-6*sqrt((mean2-mean.*mean)))
    %xlabel("Position (m)")
    %ylabel("Stress Pa")
    %legend(["mean","mean+6 standard deviations","mean-6 standard deviations"])
    %title("Nominal Spar Stress vs Position")
    %hold off
    frac=sigma/sigmaUlt;
    end
    %set up complex step
    h=1e-30;
    g=zeros([2*(N+1),N+1]);
    %calculate jacobian
    for ir=1:2*(N+1)
    comp=R;
    comp(ir)=R(ir)+complex(0,h);
    g(ir,:)= imag(subs(comp))/h;
    end
    frac=subs(R);
    %calculate outputs
    c=frac-ones(size(frac));
    ceq=[];
    jeq=[];
    j=g;
    %run convergence history
    %frac=[[],[],[],[],[],[],[]];
    %for asd=1:6
     %   ngq=asd;
     %   frac(asd,:)=subs(R);
    %end
    %calculate error
    %error=[0,0,0,0,0,0];
    %for asd=1:6
    %    error(asd)=norm(frac(asd)-frac(6));
    %end
    %plot errors
    %plot(1:6,error)
    %xlabel("Number of Gauss Quadrature Points")
    %ylabel("Norm of the error")
    %title("Norm of the error at different Gauss Quadrature points for the nominal spar")
end


