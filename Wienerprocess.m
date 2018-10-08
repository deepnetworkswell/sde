
% % % % % % % % % % % % % % % % % % 

n=1000; M=10000; % Number of time steps/trials
R=zeros(n+1,M); % Preallocate result matrix, R
for i=1:M
R(1:n+1,i)=Wiener(n); % Compute and store evaluated paths
end
expectation=sum(R,2)/M;
figure
plot(t,expectation,t,exp(9*t/8),'r-.')
title(['The mean of ' int2str(M) ' trials vs actual expectation']),...
xlabel('t'),ylabel('exp(t+0.5W(t))')
legend('Experiment','Theoretical expectation exp(9t/8)','Location','NW')



% % % % % % % % % % % % % % % % 
T=1; % Set time interval [0,T]
n=1000; % Set number of steps to compute at in [0,T]
dt=T/n; % Compute the time step
%randn('state',0) % Repeatable trials on/off
dW=sqrt(dt)*randn(1,n); % Generate array of brownian movements
W=cumsum(dW); % Sum the array cummulativley
t=0:dt:T; % Array of equal time steps
W=[0,W]; % Set first value to be zero
plot(t,W); % Plot the Wiener process
xlabel('t');ylabel('W(t)','Rotation',0);title('Sample Wiener Process')




