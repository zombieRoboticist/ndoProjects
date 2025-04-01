feasoptinocomstep=[     0.000e+00  ,       3.299e-02 ;    0.000e+00     ,     1.759e-02  ;     0.000e+00   ,      1.269e-02  ;    0.000e+00   ,      1.015e-02  ; 0.000e+00     ,   1.266e-02  ; 0.000e+00    ,     2.878e-02  ;    0.000e+00      ,  2.398e-02 ;     0.000e+00     ,    1.913e-02 ;    0.000e+00     ,   1.304e-02  ;   0.000e+00    ,     7.870e-03  ;      0.000e+00,        4.605e-03 ; 0.000e+00     ,    1.504e-03  ; 3.650e-07   ,  5.631e-05  ;   4.770e-18,       3.239e-15  ];
comstep=[  0.000e+00     ,     3.299e-02  ;       0.000e+00    ,     1.757e-02  ;      0.000e+00     ,     1.268e-02  ;      0.000e+00    ,    1.012e-02  ;     0.000e+00     ,    1.254e-02  ;     0.000e+00     ,     2.873e-02  ;      0.000e+00     ,     2.395e-02  ;        0.000e+00    ,    1.898e-02  ;       0.000e+00    ,   1.305e-02  ;     0.000e+00    ,    7.883e-03  ;     0.000e+00   ,     4.550e-03  ;    0.000e+00    ,    1.485e-03  ;     3.843e-11     ,    2.745e-08  ];
close all


%plot  feasibility graph
%plot(1:14,feasoptinocomstep(:,1))

semilogy(1:14,feasoptinocomstep(:,2))
hold on
%plot(1:13,comstep(:,1))
semilogy(1:13,comstep(:,2))
xlabel("fmincon Itterations")
title("Convergence Plot")
ylabel("First Order Optimality")
legend( "NonComplex Step First-Order Optimality","Complex Step Feasibility", "Complex Step First-Order Optimality")
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
