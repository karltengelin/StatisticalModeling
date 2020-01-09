function [th,th_cov,th_int]=inhom_poisson_est(S,t)
% inhom_poisson_est  Inhomogeneous Poisson process parameter estimation, lab2
%
% [th,th_cov,th_int]=inhom_poisson_est(S,T)
%
% S      : The process domain length.
% T      : Observed events.
% th     : Estimated model parameter values.
% th_cov : Empirical Fisher information estimate of the covariances.
% th_int : Approximate 95% confidence intervals for the parameters.
%
% Poisson process on the interval [0,S]
% Observed events: T(k), k=1,...,N(S)
% Intensity model:
%   lambda(t) = theta(1)+theta(2)*cos(2*pi*t/S)+theta(3)*sin(2*pi*t/S)

t = t(:);
n = length(t);

% Initial estimate:
th = [n/S; 0; 0];
% A few iterations of robustified Newton optimisation:
for k=1:10
  [dl,d2l] = inhom_poisson_deriv(th,S,t);
  th = adjust_estimate(th - d2l\dl,th);
end
% Fisher information based covariance estimate:
[dl,d2l] = inhom_poisson_deriv(th,S,t);
th_cov = inv(d2l);
% Approximate 95% confidence interval:
th_int = th*[1,1,1]+sqrt(diag(th_cov))*[-1.96,0,1.96];

function th=adjust_estimate(th,th0)

% The mean intensity must not be negative:
if (th(1)<0)
  th(1) = abs(th0(1))/2;
end
% The minimum intensity must not be negative:
if (th(1)^2<th(2)^2+th(3)^2)
  th(2:3) = th(2:3)/sqrt(th(2)^2+th(3)^2)*th(1)*0.99;
end
