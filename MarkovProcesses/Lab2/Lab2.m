lambda = 2;

u=rand(1,100);
y=-log(1-u)/lambda;
T=cumsum(y);
N=(1:100);
figure()
hold on
stairs(T,N)
u=rand(1,100);
y=-log(1-u)/lambda;
T=cumsum(y);
stairs(T,N)
u=rand(1,100);
y=-log(1-u)/lambda;
T=cumsum(y);
stairs(T,N)
%% 2.3.1
n=20;
lambda_est = N(n)/T(n);
alpha = norminv(0.975);
conf_lowend = lambda_est-alpha*sqrt(lambda_est/n);
conf_highend = lambda_est+alpha*sqrt(lambda_est/n);
%% 2.3.2
n=90;
lambda_est = N(n)/T(n);
alpha = norminv(0.975);
conf_lowend = lambda_est-alpha*sqrt(lambda_est/n);
conf_highend = lambda_est+alpha*sqrt(lambda_est/n);

%% 2.4
clear all
load coal.dat
T=cumsum(coal(:,5));
N=1:size(coal,1);
figure()
stairs(T,N)
lambda_est = (N(length(N))/T(length(T)))*365;
alpha = norminv(0.975);
alpha2 = norminv(0.025);
conf_lowend = lambda_est-alpha*sqrt(lambda_est/T(length(N)));
conf_highend = lambda_est+alpha*sqrt(lambda_est/T(length(N)));

% split coal data

T1 = T(100);
T2 = T(length(T))-T(101);
N1 = 100;
N2 = length(T)-100;

lambda_est1 = (N1/T1)*365;
lambda_est2 = (N2/T2)*365;

conf_lowend2 = lambda_est1 -lambda_est2+((alpha2)^2)/2*(1/N1-1/N2)-alpha*sqrt((lambda_est1/N1+lambda_est2/N2)+((alpha2)^2)/4*(1/N1-1/N2)^2);
conf_highend2 = lambda_est1 -lambda_est2+((alpha2)^2)/2*(1/N1-1/N2)+alpha*sqrt((lambda_est1/N1+lambda_est2/N2)+((alpha2)^2)/4*(1/N1-1/N2)^2);

figure()
stairs(coal(:,3)+coal(:,4)/365,N)
figure()
stairs(coal(:,3)+coal(:,4)/365,cumsum(coal(:,6)))

%% 2.5
clear all
S = 0.5;
lambda = 500;

y= porand(lambda*S);

i=1;
counter = 1;

while counter < y
    K = rand(1,1);
    T = rand(1,1);
    if K+T < 1
        xtemp(counter) = K;
        ytemp(counter) = T;
        counter = counter+1;
    end
end
figure()
plot(xtemp,ytemp,'x')






% lambda = 0.5;
% 
% Y = 1:100;
% X = zeros(1,length(Y));
% 
% for i = 1:length(X)
%     
% X(1,i) = porand(lambda)+1;
% 
% end
% 
% X = cumsum(X);
% figure()
% stairs(X,Y)








%% 3.1
clear all
u=0:0.01:10;
S = 10;
gamma = 8;

y= porand(gamma*S);
T = sort(rand(y,1)*S);

figure()
plot(u,4*(1-cos(u.^2/4)),T,0*T,'x')

%% 3.2.1
S = 50; 
theta = 0.5*[20, 15, 0];
t = inhom_poisson_simulate(theta,S);
u = 0:S/(length(t)-1):S;
figure()
plot(u,theta(1,1) + theta(1,2)*cos(2*pi*t/S) + theta(1,3)*sin(2*pi*t/S),t,0*t,'x')
%% 3.2.2
S = 10; 
theta = [20, 15, 0];
t = inhom_poisson_simulate(theta,S);
[tau, tau2, tau3] = inhom_poisson_est(S,t);
u = 0:S/(length(t)-1):S;
figure()
plot(u,theta(1,1) + theta(1,2)*cos(2*pi*t/S) + theta(1,3)*sin(2*pi*t/S),t,0*t,'x')


%% 3.2.3
S = 10; 
theta = 0.5*[20, 15, 0];
t = inhom_poisson_simulate(theta,S);
[tau, tau2, tau3] = inhom_poisson_est(S,t);
u = 0:S/(length(t)-1):S;
figure()
plot(u,theta(1,1) + theta(1,2)*cos(2*pi*t/S) + theta(1,3)*sin(2*pi*t/S),t,0*t,'x')