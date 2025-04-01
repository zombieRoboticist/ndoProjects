function [c] = suragate(x,z)

% test the gp library
%clear all;
%close all;

% this is needed so Matlab knows where to find the library if the gpml
% directory is not the working directory; set mydir to the location of
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

% define the function to be sampled
fe = @(x) -1*obj(x);

% sample the function at mid point and bounds
y=zeros([length(x),1]);
for k=1:length(x)
    y(k)=fe(x(k,:));
end

% set the squared exponential covariance function
covfunc = @covSEiso; % {@covMaterniso, 1}; % 
hyp.cov = [log(1); log(10000)]; % first component is log(l) and second is log(sigma)

% set the likelihood function to Gaussian
likfunc = @likGauss;
sn = 1000; %1e-16; % this is the noise level
hyp.lik = log(sn);

% maximize the likelihood function to find the hyperparameters
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, x, y);

%fz=@(zt) gp(hyp, @infExact, [], covfunc, likfunc, x, y, z);
% now sample the GP for plotting purposes.  In project 3, the variable z
% will be a single design provided to the objective function by fmincon.
% Here, we want to plot the GP surrogate, so we need to sample many z
% values.
%z = [linspace(0.1, 1.5, 125)',linspace(0, .3, 125)',linspace(3*pi/30, 8*pi/30, 125)'];
m = gp(hyp, @infExact, [], covfunc, likfunc, x, y, z');
c=m;
%zt=1:125;
% ...and plot
%f = [m+1.96*sqrt(s2); flipdim(m-1.96*sqrt(s2),1)];
%fill([zt'; flipdim(zt',1)], f, [7 7 7]/8)
%hold on;
%plot(zt, m);
%plot(x, y, 'ks','MarkerFaceColor','k');
%plot(zt, fe(z), 'k--')
%axis([0 7.5 -3000 3000]);
%hold off;
