%------Firstly we initialize all the matrices and vectors----%
Alpha = [0 2/5 1/5 0 0; 0 0 3/4 1/4 0; 1/2 0 0 1/2 0; 0 0 1/3 0 2/3; 0 1/3 0 1/3 0];
w = Alpha*ones(size(Alpha,2),1);
r = max(w);
D = diag(w);
P_targets = D\Alpha;

% number of iterations or expected running time:
Tmax = 20; 

% Transition probability matrix:

P_transprob = zeros(size(P_targets));
for i = 1:length(P_transprob)
    for j = 1:length(P_transprob)
        if i ~= j
            P_transprob(i,j) = Alpha(i,j)/max(w);
        end
    end
    P_transprob(i,i) = 1-sum(P_transprob(i,:));
end

%% -- From here we specify the simulation -- %%

nbr_of_part = 100;                      % number of particles
nbr_of_sim = 10000;                     % number of simulations
times = zeros(1,nbr_of_part);           % initializing the vector in which we store times
first_return = zeros(1,nbr_of_sim);
start = 2;                              % We specify our starting node: 1-o, 2-a, 3-b, 4-c, 5-d
finish = 2;                             % We specify our finishing node: 1-o, 2-a, 3-b, 4-c, 5-d

for iii = 1:nbr_of_sim
    for ii = 1:nbr_of_part
    %---- Simulate the particle moving around ----%

    x = zeros(5,Tmax);                  % We specify the network in matrix form
    x(start,1) = 1;                     % setting starting node = 1        
    n = start;                          % pre determining our initial node               
    time = zeros(1,Tmax-1);             % initialize a time vector in which to store times
    %%
        for i = 2:Tmax
            random = rand;              % Here we randomize where the particle jumps to next

                if random <= sum(P_transprob(n,1:1))
                    n = 1;
                elseif random <= sum(P_transprob(n,1:2))
                    n = 2;
                elseif random <= sum(P_transprob(n,1:3))
                    n = 3;    
                elseif random <= sum(P_transprob(n,1:4))
                    n = 4;
                elseif random <= sum(P_transprob(n,1:5))
                    n = 5;
                end
            x(n, i) = 1;                % updating the state vector 
            time(i-1) = -log(rand)/r;   % assigning a time to the taken step

            if x(finish,i) == x(start,1) && x(start,i-1) ~= x(start,1)  %if we have reached our finish node given that we have left our starting node we break the loop so that no more time is added
                break
            end
            
        end
        times(ii) = sum(time);          % storing the total time it took for the particle from start to finish
    end
    first_return(iii) = min(times);     % taking the smallest time it took for a particle to get from start to finish
end
mean(first_return)                      % taking the mean time of all the first returned particles


