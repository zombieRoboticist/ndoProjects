maxItter =200;
changeTol=10^-9;
accuracyTol=5*10^-2;

% Initialize x
%r20 = [.1,.8,1.5];
%a0=[0,.15,.3];
%w0=[3*pi/30,5*pi/30,8*pi/30];
%x=zeros([length(r20)*length(a0)*length(w0),3]);
%for i=1:length(r20)
%    for j=1:length(a0)
%        for k=1:length(w0)
%            x(length(a0)*length(w0)*(i-1)+length(w0)*(j-1)+k,:)=[r20(i),a0(j),w0(k)];
%        end
%    end
%end
x=lhsdesign(200,3);
x(:,1)=.1+x(:,1)*1.4;
x(:,2)=0+x(:,2)*.3;
x(:,3)=3*pi/30+x(:,3)*5*pi/30;
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
fe = @(x) -1*obj(x,800);

% set the squared exponential covariance function
covfunc =  {@covMaterniso, 5}; % 
hyp.cov = [log(.4); log(.75)]; % first component is log(l) and second is log(sigma)

% set the likelihood function to Gaussian
likfunc = @likGauss;
sn = .4; %1e-16; % this is the noise level
hyp.lik = log(sn);

%set fmincon options
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
lb=[.1,0,3*pi/30];
ub=[1.5,.3,8*pi/30];

a0=[(ub(1)+lb(1)-.5)/2; (ub(2)+lb(2)+.1)/2;(ub(3)+lb(3)+.45)/2];
count=0;
% sample the function at xs
    y=zeros([length(x),1]);
    for k=1:length(x)
        y(k)=fe(x(k,:));
    end
errorobj=[];
erroritter=[];
while count<maxItter


    
    
    % maximize the likelihood function to find the hyperparameters
    hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, x, y);


    fx=@(z) gp(hyp, @infExact, [], covfunc, likfunc, x, y, z');
    %run optimization function
    a=fmincon(fx,a0,[],[],[],[],lb,ub,[],options);
    %convergend stopping conditions
    if  abs(fx(a)+fe(a))<accuracyTol
        if all(abs(fx(a)-fx(a0))<changeTol)
            break
        end
    end
    %update variables
    errorobj=[errorobj;fx(a)+fe(a)];
    erroritter=[erroritter;fx(a)-fx(a0)];
    x=[x;a'];
    y=[y;fe(a)];
    a0=a;
    count=count+1;
end
a