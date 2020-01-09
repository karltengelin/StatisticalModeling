function [dl,d2l]=inhom_poisson_deriv(th,S,t)
% inhom_poisson_deriv  Inhomogeneous Poisson process likelihood derivatives
%
% [dl,d2l]=inhom_poisson_deriv(th,S,T)
%
% th     : Model parameter values.
% S      : The process domain length.
% T      : Observed events.
% dl     : The derivatives of the negative log-likelihood.
% d2l    : The 2nd order derivatives of the negative log-likelihood.
%
% Poisson process on the interval [0,S]
% Observed events: T(k), k=1,...,N(S)
% Intensity model:
%   lambda(t) = theta(1)+theta(2)*cos(2*pi*t/S)+theta(3)*sin(2*pi*t/S)
% Negative log-likelihood:
%   l(th) = Lambda(S) + log(N(S)!) - sum(log(lambda(T)))

t = t(:);
B = [ones(length(t),1), cos(2*pi*t/S), sin(2*pi*t/S)];
lambda = inhom_poisson_lambda(th,S,t);
dl = [S;0;0]-B'*(1./lambda);
d2l = B'*(((1./lambda.^2)*[1,1,1]).*B);
