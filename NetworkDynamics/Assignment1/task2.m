clear all
close all
load -ASCII twitter.mat
W = spconvert(twitter);
index = find(~sum(W,2)); % find any out degress that sum to zero (i.e sinks) and return index
W(index,index) = 1; %implement a self loop in the sink so we don't lose information
beta = 0.15;

sz = max(size(W)); 
W(sz,sz) = 0; % Since W is a sparse matrix we need to force matlab to return a square matrix

D_inv = diag(sum(W,2).^(-1)); % Here we create the D-inverse in order to calculate the P matrix in the next step
P = D_inv*W;
%% a
mu = ones(length(W),1); %Creating the mu-vector
mu_i = mu.^(-1); %Creating the inverse of the mu-vector

z_old = 0; %initializing z_old, z and P0
z = beta*mu;
P0 = beta*mu;

% This while loop updates the z according to equation 2.5 in lecture notes until the 
% difference between the updated value and the last value is small (i.e until convergence)
while norm(z-z_old,Inf) > 10^(-6) 
z_old = z;
P0 = (1-beta)*P'*P0;
z = z + P0;
end

%now we extract the 5 largest and 5 smalles values from the converged z-vector
[values_central central_nodes] = maxk(z,5);
[values_noncentral noncentral_nodes] = mink(z,5);
%displaying most central nodes and least central nodes
central_nodes
noncentral_nodes

%% b and c

%First we choose stubborn nodes
stubb = [9 214];

% check which nodes that are regular given the initialized stubborn nodes
reg = setdiff(1:length(P),stubb);

%initialize the startvalue for all non stubborn nodes
commonstartval = 0.5;

%creating Q and E-matrix
Q = P(reg,reg);
E = P(reg,stubb);

%setting the value of the stubborn nodes
u_1 = [1;0];
%setting the start values of the regular nodes
x_1 = commonstartval*ones(length(Q),1);
% choosing 10 random nodes to observe
random = randperm(length(P),10);


iter = 250; % number of iterations
plotmat = zeros(length(Q),iter);    % creating the matrix where the value in in each node is stored for each
                                    % iteration, will later be used to plot the value of each random node over
                                    % time
k = 0;
while k < iter %looping over number of iterations
    pause(0.01);
    x_1 = Q*x_1+E*u_1; %updating the values of all the regular nodes
    k = k+1;
    x = [u_1;x_1];
    %bar(random,x(random))
    hist(x)
    plotmat(1:stubb(1)-1,k) = x_1(1:stubb(1)-1); %storing the value of each regular node for every iteration
    plotmat(stubb(1)+1:stubb(2)-1,k) = x_1(stubb(1)+1:stubb(2)-1);
    plotmat(stubb(2)+1:end,k) = x_1(stubb(2)+1:end);
end

%setting the startvalues in the plot-matrix
plotmat(:,1) = commonstartval*ones(length(plotmat),1);
plotmat(stubb(1),:) = u_1(1)*ones(1,iter);
plotmat(stubb(2),:) = u_1(2)*ones(1,iter);

% Plotting
figure
hold on
%stubborn nodes
plot(plotmat(stubb(1),:))
plot(plotmat(stubb(2),:))
%non-stubborn nodes
for i = 1:length(random)
plot(plotmat(random(i),:))
end
hold off
