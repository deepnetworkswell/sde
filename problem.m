function problem
clear all
clc
field1 = 'TimeInterval';  value1 = 10;
field2 = 'Steps';         value2 = 2^17;   % Number of steps
field3 = 'X0';            value3 = 10;     % Initial condition
field4 = 'V0';            value4 = 0;
m=1 ;    gamma=10 ;    k=10 ;     Vs=10 ; 
omega2=k/m ; alpha=k*Vs/m ; nu=gamma/m;  beta=1/m ;
Constants= struct('m',m,'gamma',gamma,'k',k,'Vs',Vs,'omega2',omega2,'alpha',alpha,'nu',nu,'beta',beta,'noiseamp',10);
Params = struct(field1,value1,field2,value2,field3,value3,field4,value4,'Constants',Constants);
T=Params.TimeInterval;
steps=Params.Steps;
dt=T/steps; % Set stepsize
nruns=3;
steps=steps+1;
X=zeros(steps,nruns);V=zeros(steps,nruns);
for k=1:nruns
    [X(:,k),V(:,k)]=onerun(Params);
end
expectationX=sum(X,2)/nruns;
expectationV=sum(V,2)/nruns;
springStretch_simulation=expectationX(end)-Vs*T  %spring stretch in last step
springStretch_theory=2*nu*Vs/omega2
figure
plot(0:dt:T,expectationX,'-',0:dt:T,expectationV,'g.');
title('Average Numerical solution of an SDE using the Euler-Maruyama Method')
xlabel('t');
legend('X','V');
end
%%

function [X,V]=onerun(Params)
T=Params.TimeInterval;
steps=Params.Steps;
X0=Params.X0;
V0=Params.V0;
dt=T/steps; % Set stepsize
dW=sqrt(dt)*randn(1,steps); % Array of Brownian movements
W=cumsum(dW); % Wiener process
Y1=zeros(steps,1);Y2=zeros(steps,1); % Create array for results
Y1old=X0; % Initialise Xold
Y2old=V0; % Initialise Xold
% Set the Parameters*************
m=Params.Constants.m ;    gamma=Params.Constants.gamma ;    k=Params.Constants.k ;    
Vs=Params.Constants.Vs ; omega2=Params.Constants.omega2 ; alpha=Params.Constants.alpha ; 
nu=Params.Constants.nu;  beta=Params.Constants.beta ; noiseamp=Params.Constants.noiseamp; 
%********************************
for i=1:steps
    Y1(i) = Y1old + Y2old * dt;
    Y2(i) = Y2old - nu*Y2old*dt - omega2*Y1old*dt ...
        + alpha/2*(2*i+1)*dt*dt+ noiseamp*beta*dW(i)+ 000* heaviside(Y2old)*dt;    %%alpha*dt*dt
    Y1old=Y1(i); % Update Xold
    Y2old=Y2(i); % Update Vold
end
Y1=[X0;Y1];Y2=[V0;Y2];
% figure;
% plot(0:dt:T,Y1,'-',0:dt:T,Y2,'r.');
% title('Numerical solution of an SDE using the Euler-Maruyama Method')
% xlabel('t');
% legend('X','V');
X=Y1;
V=Y2;
end
