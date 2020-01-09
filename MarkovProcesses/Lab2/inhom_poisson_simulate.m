function t=inhom_poisson_simulate(th,S)
% inhom_poisson_simulate  Simulate an inhomogeneous Poisson process
%
% T=inhom_poisson_simulate(th,S)
%
% th     : The model parameter values.
% S      : The process domain length.
% T      : Simulated event time points.
%
% Poisson process on the interval [0,S]
% Events: T(k), k=1,...,N(S)
% Intensity model:
%   lambda(t) = theta(1)+theta(2)*cos(2*pi*t/S)+theta(3)*sin(2*pi*t/S)

lambda_sup = th(1)+sqrt(th(2)^2+th(3)^2);
n = porand(lambda_sup*S);
t = sort(rand(n,1)*S);
[lambda,Lambda] = inhom_poisson_lambda(th,S,t);
t = t(rand(n,1)<=lambda/lambda_sup);
