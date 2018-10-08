function [expW] = Wiener(n,T)
if nargin==1, T=1; end % Sets T=1 if only 1 argument is input
dt=T/n;
dW=sqrt(dt)*randn(1,n);
W=cumsum(dW); % Generates Brownian path
t=[0:dt:T]; W=[0,W];
expW=exp(t+0.5*W); % Evaluates the function along the path
end