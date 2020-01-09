%% -- initializing the neccessary matricex and vectors -- %%
Alpha = [0 2/5 1/5 0 0; 0 0 3/4 1/4 0; 1/2 0 0 1/2 0; 0 0 1/3 0 2/3; 0 1/3 0 1/3 0];
w = Alpha*ones(size(Alpha,2),1);
r = max(w);
D = diag(w);
P_targets = D\Alpha;

P_transprob = zeros(size(P_targets));
for i = 1:length(P_transprob)
    for j = 1:length(P_transprob)
        if i ~= j
            P_transprob(i,j) = Alpha(i,j)/max(w);
        end
    end
    P_transprob(i,i) = 1-sum(P_transprob(i,:));
end

%% 1b
[V,D] = eig(P_transprob');              % Calculating values and associated eigenvectors for transition probability matrix
lambda = diag(D);                       % Putting al eigenvalues in a diagonal matrix
[~,i] = sort(lambda,1,'descend');       % Sorting eigenvalues to find the biggest
pi_hat = V(:,i(1));                     % setting pi_hat as the eigenvector with largest eigenvalue
pi_hat = pi_hat/sum(pi_hat);            %normalizing pi_hat
theoretic_return_node_a = 1/pi_hat(2)   % Returning the theoretical return time, we set pi_hat(2) since node a corresponds to node number 2

%% 1d
S = [5];                                        % Defining our subset S containing the nodes we aim to hit
R = setdiff(1:5,S);                             % Removing rows and columns associated with S in trans. prob. matrix 
P_new = P_transprob(R,R);
x = (eye(length(R))-P_new)\ones(length(R),1);   % Doing the matrix calculations to solve the linear equation system, here x(i) will contain the time it takes to reach S from i 
hitting_time_node_o = x(1)                      % Returning the hitting time for node o to node d

