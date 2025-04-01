    %%run fmincon first
    % where you saved the gpml directory
mydir = '../gpml-matlab-v36/';
addpath(mydir(1:end-1))
addpath([mydir,'cov'])
addpath([mydir,'doc'])
addpath([mydir,'inf'])
addpath([mydir,'lik'])
addpath([mydir,'mean'])
addpath([mydir,'prior'])
addpath([mydir,'util'])
%define actual generalized model
    fx=@(z) gp(hyp, @infExact, [], covfunc, likfunc, x, y, z);

xx=3:6/100:9;
yx=zeros(size(xx));
ys=zeros(size(xx));
rx=zeros(size(xx));
rs=zeros(size(xx));
%evaluate objective and generalized models
for i=1:100
    yx(i)=obj([.8,.058,xx(i)*pi/30],800);
    ys(i)=-1*fx([.8,.058,xx(i)*pi/30]);
    rx(i)=obj([a(1),a(2),xx(i)*pi/30],800);
    rs(i)=-1*fx([a(1),a(2),xx(i)*pi/30]);
end
%plot results
plot(xx,yx)
hold on 
%plot(xx,ys)
%legend(["Objective Function", "Suragate Model"])
%title("Objective function and Optimized Suragate Model at r2=.8m and a1=.058 rad")
title("Objective Function at r2=.8m and a1=.058 rad")
ylabel("Standard Deviation of the Dynamics Model")
xlabel("Speed (rotations per second)")
hold off
pause
plot(xx,rx)
hold on
plot(xx,rs)
legend(["Objective Function", "Suragate Model"])
title("Objective function and Suragate Model on Optimized Slice")
ylabel("Standard Deviation of the Dynamics Model")
xlabel("Speed (rotations per second)")
hold off

%%
T=10:1000;
sigma=zeros(size(T));
for i=1:numel(T)
    sigma(i)=obj([.8,.058,6.5*pi/30],T(i));
end
plot(T,sigma)
hold on
title("Convergence of Objective Function over Movement Cycles")
ylabel("Standard Deviation of the Dynamics Model")
xlabel("Number of Hills")
hold off

