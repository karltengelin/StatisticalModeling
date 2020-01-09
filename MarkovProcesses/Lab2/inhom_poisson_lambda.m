function [lambda,Lambda]=inhom_poisson_lambda(th,S,t)
% inhom_poisson_lambda  Inhomogeneous Poisson process intensity
%
% [lambda, Lambda]=inhom_poisson_lambda(th,S,t)
%
% th     : The model parameters.
% S      : Process domain interval length.
% t      : Time points.
% lambda : The intensity at the given time points.
% Lambda : The integrated intensity.
%
% Poisson process on the interval [0,S]
% Events: T(k), k=1,...,n
% Intensity model:
%   lambda(t) = theta(1)+theta(2)*cos(2*pi*t/S)+theta(3)*sin(2*pi*t/S)

lambda = th(1)+th(2)*cos(2*pi*t/S)+th(3)*sin(2*pi*t/S);
Lambda = th(1)*S;
