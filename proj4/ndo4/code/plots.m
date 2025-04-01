feasoptinocomstep=[    0.000e+00    ,     4.712e-02 ;      0.000e+00 ,     3.193e-02  ;    0.000e+00     ,     2.027e-02  ;  0.000e+00    ,    1.155e-02  ;   0.000e+00    ,     9.005e-03  ;  0.000e+00     ,     1.209e-02  ;   0.000e+00     ,     1.461e-02  ;   0.000e+00     ,   1.866e-02  ;  0.000e+00     ,     2.231e-02  ;   0.000e+00     ,    1.833e-02  ;  0.000e+00     ,    1.372e-02  ;  0.000e+00     ,    7.524e-03  ;   0.000e+00    ,    3.030e-03  ; 1.532e-06     ,    1.337e-04  ; 6.503e-11     ,  2.741e-07  ];

comstep=[  0.000e+00     ,     3.299e-02  ;       0.000e+00    ,     1.757e-02  ;      0.000e+00     ,     1.268e-02  ;      0.000e+00    ,    1.012e-02  ;     0.000e+00     ,    1.254e-02  ;     0.000e+00     ,     2.873e-02  ;      0.000e+00     ,     2.395e-02  ;        0.000e+00    ,    1.898e-02  ;       0.000e+00    ,   1.305e-02  ;     0.000e+00    ,    7.883e-03  ;     0.000e+00   ,     4.550e-03  ;    0.000e+00    ,    1.485e-03  ;     3.843e-11     ,    2.745e-08  ];
%

%plot  feasibility graph
%plot(1:14,feasoptinocomstep(:,1))

semilogy(0:14,feasoptinocomstep(:,2))
hold on
%plot(1:13,comstep(:,1))
%semilogy(1:13,comstep(:,2))
xlabel("fmincon Itterations")
title("Convergence Plot")
ylabel("First Order Optimality")
xlabel("Number of Iterations")
legend( "First-Order Optimality")
hold off

%%
% calculate convergence study

%initiialize the design constants
L=7.5;
m=500/2;
%N=50;
rho=1600;
E=70000000000;
sigmaUlt=600000000;
start=10;
stop=200;

x=start:stop;
q=zeros(size(x));
z=1;
for i=start:stop
    R=putR(.015*ones([i+1,1]),.05*ones([i+1,1]));
    [c,j,a,b]=getC(R,L,m,i,E,sigmaUlt);
    q(z)=c(i+1)+1;
    z=z+1;
end

%%
%plot convergence study

semilogy(x,q)
xlabel("Number of Elements")
ylabel("log(normalized stress)")
title("Normalized Stress at Tip of Spar")
