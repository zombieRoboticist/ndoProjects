%initialize the design constants
L=7.5;
m=500/2;
N=80;
rho=1600;
E=70000000000;
sigmaUlt=600000000;

R=putR(.01*ones([N+1,1]),.05*ones([N+1,1]));
%[rin,rout]=getR(R);
%initialize function
vol=@(R) calcVol(R,L,N);
%initialize c nonlcon
chat=@(R) getC(R,L,m,N,E,sigmaUlt,2);

%set fmincon options
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);

%initialize lb,ub, a,b

lb=putR(.01*ones([N+1,1]),.0125*ones([N+1,1]));
ub=putR(.0475*ones([N+1,1]),.05*ones([N+1,1]));
bineq=-.0025*ones([N+1,1]);
Aineq=zeros([N+1,2*(N+1)]);
for i=1:length(bineq)
    Aineq(i,2*i-1)=1;
    Aineq(i,2*i)=-1;
end
%run optimization function
a=fmincon(vol,R,Aineq,bineq,[],[],lb,ub,chat,options);
%plot output design function
[ain,aout]=getR(a);
x=0:L/N:L;
plot(x,ain,'b')
hold on
plot(x,aout,'r')
plot(x,-1*aout,'r')
plot(x,-1*ain,'b')
plot(x,zeros(size(x)),"-k")
xlabel("Length (m)")
ylabel("Radius (m)")
legend('Inner Radius','Outer Radius')
title("Optimized Spar Design")
hold off


