%%This file runs the optimizer and calc flux functions and outputs graphs
%%of the resulting outputs.

%% get single optimiztion 

%define Nx, Ny for the optimizer mesh
Nx=500;
Ny=Nx;
%define initial conditions
a0=[.012;ones(15,1)*.001];
%run optimization code
a=run_opt(Nx,Ny,a0);

%% plot function (Note: Nx and vector a of the coefficients used to define 
% the height function must de defined prior to running this code. This can 
% be done by runing the get single optimization section)
%define vector of positions
x = [0:L/Nx:L].';
%plot the height function as a function of x
plot(x,calcHeight(a,.05,x),'red')
hold on
%plot(x,calcHeight(a0,L,x), 'blue')
%add y axis label
ylabel("Height (m)")
%add x axis label
xlabel("Position (m)")
%add plot title
title("Sectional View of the Extruded Heat Exchanger")
hold off


%[flux,T,dTdx,xy] = CalcFlux(.05, calcHeight(a,.05,x), Nx, Ny, 20, 20, 90);
%flux

%% plot Convergence (Note: run get single optimization code and put the 
% correct terminal output into feas)
%define array of itteration numbers
itter=0:1:13;
%define matrix of feasablity (first col) and first order optimality (second
%col) using printed output from fmincon
feas=[ 9.223e-03    ,     2.036e-02  ;   3.469e-18    ,     1.710e-02  ;   1.041e-17     ,     5.476e-03  ;    5.204e-18     ,    1.431e-03  ;    1.214e-17     ,     3.199e-04  ;    8.674e-18    ,   3.389e-05  ;     2.082e-17     ,   3.985e-05 ;    1.388e-17     ,   1.634e-05  ;   3.469e-17     ,     1.800e-05  ;   4.857e-17     ,     4.567e-05  ;    3.469e-17    ,    1.536e-05  ;    3.469e-17    ,     3.119e-05  ;   2.776e-17    ,   1.394e-06  ;   2.776e-17     ,    2.575e-19   ];
%plot log(feasibility) vs itteration number 
semilogy(itter,feas(:,1),'blue')
hold on
%plot log(first order optimality) vs itteration number
semilogy(itter,feas(:,2),'red')
%add chart title
title("Feasibility and Convergance History")
%label x axis
xlabel("Number of Itterations")
%label y axis
ylabel("log(Feasibility, First Order Optimality)")
%add legend
legend(["Feasibility","First-Order Optimality"])
hold off

%% run convergence study (Note: vector a of the coefficients used to define
% the height function must de defined prior to running this code)
%setup initial and final nxs
nstart=10;
nfinal=1000;
%create vector of nxs to try
nx=[nstart:nstart:100,200:100:nfinal];
%initialize vector to store fluxes
fluxx=[];
for i=1:max(size(nx))
    %calculate x for the current nx
    x = [0:L/nx(i):L].';
    %calculate the heat flux for the current nx
    [flux,T,dTdx,xy] = CalcFlux(.05, calcHeight(a,.05,x), nx(i), nx(i), 20, 20, 90);
    %append the current heat flux to the end of the vector storing the heat
    %fluxes
    fluxx=[fluxx,flux];
end


%% plot convergence study (Note: run convergence study code before running 
% this part)
%plot flux vs nx
plot(nx,fluxx)
hold on
%plot lower bound of delta
plot(nx,ones(size(fluxx))*(fluxx(max(size(fluxx)))*.95),"--b")
%plot upperbound of delta
plot(nx,ones(size(fluxx))*(fluxx(max(size(fluxx)))*1.05), "--b")
%add title
title("Convergence Study for the Optimized Heat Exchanger Design")
%add x axis label
xlabel("Number of Segments (Nx =Ny)")
%add y axis labelflux
ylabel("Flux/Unit Length (W/m)")
%add legend
legend(["Convergence plot","flux( Nx=1000)+/- 5%"])
hold off