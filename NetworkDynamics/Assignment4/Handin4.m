%% Preliminary parts
% Epidemic on a known graph
n = 500;    % Number of nodes
W = zeros(n);
W = W + diag(ones(n-1,1),1); % add ones on the +1 off-diagonal
W = W + diag(ones(n-1,1),-1); % add ones on the -1 off-diagonal
W = W + diag(ones(n-2,1),2); % add ones on the +2 off-diagonal
W = W + diag(ones(n-2,1),-2); % add ones on the -2 off-diagonal
W = W + diag(ones(1,1),n-1); % add ones on the n-1 off-diagonal
W = W + diag(ones(1,1),1-n); % add ones on the -n+1 off-diagonal
W = W + diag(ones(2,1),n-2); % add ones on the n-2 off-diagonal
W = W + diag(ones(2,1),2-n); % add ones on the -n+2 off-diagonal
G = sparse(W); % transform it into a sparse matrix
%%
X0 = zeros(n,1);    % Node states
S = 0;  % Susceptible state
I = 1;  % Infected state
R = 2;  % Recovered state
V = 3;  % Vaccinated state
beta = 0.3; % Infection probability
rho = 0.7;  % Recovery probability
%% Randomly choose 10 infected patient zero
patientszero = randi(n,10,1);
X0(patientszero) = 1;
%%
N = 100;    % Number of trials.
avers = zeros(16,3);    % Placeholder for the average number of nodes in each state.
avers(1,:) = N*[500-10 10 0];   % Fill up first row with the initial state.
newinfected = zeros(16,1);  % Placeholder for the number of newly infected nodes.
newinfected(1) = N*10;  % Fill first elements with initial condition.
for i = 1:N     % We'll do 100 trails.
    X = X0; % Reset the state distribution to the initial state.
    for t = 1:15    % t = 0 is our initial state, we iterate from there.
        infect = 0; % Reset the number of infected this time step.
        Xtmp = (X == I);    % Find all infected nodes.
        m = G*Xtmp; % Find each node's number of infected neighbors.
        for k = 1:n     % Loop through all nodes.
            u = rand;   % To simulate random events.
            switch X(k) % Depending on the current state of node k, we'll do different things.
                case 0  % For a susceptible node, there's a chance of infection dependent on its number of infected neighbors.
                    if(u < (1-(1-beta)^(m(k))))
                    X(k) = 1;   % Node k becomes infected.
                    infect = infect + 1;    % Another node has been infected this time step.
                    end
                case 1  % An infected node might recover.
                    X(k) = (u < rho) + 1;
            end
        end
        newinfected(t+1) = newinfected(t+1) + infect;   % Store this trials results in the corrensponding time step.
        suscept = sum(X == S);
        infect = sum(X == I);
        recover = sum(X == R);
        avers(t+1,:) = avers(t+1,:) + [suscept infect recover]; % Store the number of nodes in each state.
    end
end

% Calculate the averages over N trials.
avers = avers/N;
newinfected = newinfected/N;    

% Plot the results.
plot(0:15,avers);
figure;
plot(0:15,newinfected);
%% 1.2  Generate a random graph
k = 14; % Choose an average degree.
c0 = k/2;   % Number of links to be added each time step.
k0 = k+1;   % Number of nodes in the initial complete graph.
if(mod(k0,1) < .5)  % If k is not an even integer, we need to round k0;
    k0 = floor(k0);
else
    k0 = ceil(k0);
end
W = ones(k0);   % Create an initial complete graph.
N = 950;    % Choose the final number of nodes.
% val = 1;
verts = size(W,1);  % How many nodes we currently have.

for t = 1:N-k0
    if(mod(c0,1) ~= 0)  % If k is not an even integer, we need to choose c for this iteration.
        u = rand;
        if(u < mod(c0,1))   % If the remainder is big, we want to round up more often than not.
            c = ceil(c0);
        else    % If the remainder is small, we want to round down more often.
            c = floor(c0);
        end
    else% If k is an even integer, simply use c0.
        c = c0;
    end
    W = [W zeros(size(W,2),1);  % Increase the size of W in preparation of a new node.
        zeros(1,size(W,2)+1)];
    
    w = sum(W,2);   % Node out-degrees
    P = w/sum(w);   % Node preference probabilities.
    for j = 1:c
        neighbor = randsample(1:verts+1, 1, true, full(P)); % Choose an integer dependent on the probability distribution P
        P(neighbor) = 0; % Reset the P of a chosen node, cant have more than one node to the new node.
        W(verts+1,neighbor) = 1;    % Create undirected link.
        W(neighbor,verts+1) = 1;
    end
    verts = verts + 1;  % Increase the number of nodes currently in graph.
end
check = mean(W*ones(size(W,2),1));  % Check average degree for the finished graph.
%%


%% 2 Simulate a pandemic without vaccination
k = 6; % Choose average degree.
c0 = k/2;
k0 = k+1;
if(mod(k0,1) < .5)
    k0 = floor(k0);
else
    k0 = ceil(k0);
end
W = ones(k0);
n = 500;    % Choose number of nodes of finished graph
verts = size(W,1);

for t = 1:n-k0
    if(mod(c0,1) ~= 0)
        u = rand;
        if(u < mod(c0,1))
            c = ceil(c0);
        else
            c = floor(c0);
        end
    else
        c = c0;
    end
    W = [W zeros(size(W,2),1);
        zeros(1,size(W,2)+1)];
    
    w = sum(W,2);
    P = w/sum(w);
    for j = 1:c
        neighbor = randsample(1:verts+1, 1, true, full(P));
        P(neighbor) = 0;
        W(verts+1,neighbor) = 1;
        W(neighbor,verts+1) = 1;
    end
    verts = verts + 1;
end
check = mean(W*ones(size(W,2),1));
%%  Repeated for easy access. Can be ignored if one's going from top to bottom.
X0 = zeros(n,1);
S = 0;
I = 1;
R = 2;
V = 3;
beta = 0.3;
rho = 0.7;

patientszero = randi(n,10,1);
X0(patientszero) = 1;
%%
G = sparse(W);  % Make the simulation quicker by making W sparse.
N = 100;
avers = zeros(16,3);
avers(1,:) = N*[500-10 10 0]; % All but ten infected people are initially susceptible.
newinfected = zeros(16,1);
newinfected(1) = N*10;
for i = 1:N
    X = X0;
    for t = 1:15
        infect = 0;
        Xtmp = (X == I);
        m = G*Xtmp;
        for k = 1:n
            u = rand;
            switch X(k)
                case 0
                    if(u < (1-(1-beta)^(m(k))))
                    X(k) = 1;
                    infect = infect + 1;
                    end
                case 1
                    X(k) = (u < rho) + 1;
            end
        end
        newinfected(t+1) = newinfected(t+1) + infect;
        suscept = sum(X == S);
        infect = sum(X == I);
        recover = sum(X == R);
        avers(t+1,:) = avers(t+1,:) + [suscept infect recover];
    end
end

avers = avers/N;
newinfected = newinfected/N;

plot(0:15,avers);
figure;
plot(0:15,newinfected);
%%


%% 3 Simulate a pandemic with vaccination
Vacc = [0 5 15 25 35 45 55 60 60 60 60 60 60 60 60 60]; % How large part of the population should be vaccinated by step t.
N = 100;
G = sparse(W);
avers = zeros(16,4);    % Make room for averaging number of vaccninated.
avers(1,:) = N*[500-10 10 0 0];
newinfected = zeros(16,1);
newvaccinated = zeros(16,1);
newinfected(1) = N*10;
for i = 1:N
    X = X0;
    for t = 1:15
        infect = 0;
        if(Vacc(t+1)-Vacc(t) > 0)   % If there's need to vaccinated before next time step.
            tovacc = round(Vacc(t+1)*n/100) - sum(X == V); % How many are to be vaccninated.
            notvaccinated = find(X ~= V); % Find all non-vaccinated nodes.
            bigness = length(notvaccinated);
            inds = notvaccinated(randi(bigness,tovacc,1));  % Find the indices of those who will be vaccinated.
            X(inds) = V;    % Vaccinate people.
        else
            tovacc = 0; % For later plotmaking.
        end
        
        Xtmp = (X == I);
        m = G*Xtmp;
        for k = 1:n
            u = rand;
            switch X(k)
                case 0
                    if(u < (1-(1-beta)^(m(k))))
                    X(k) = 1;
                    infect = infect + 1;
                    end
                case 1
                    X(k) = (u < rho) + 1;
            end
        end
        
        newinfected(t+1) = newinfected(t+1) + infect;
        newvaccinated(t+1) = newvaccinated(t+1) + tovacc;   % Store the vaccination rate for later plotmaking.
        suscept = sum(X == S);
        infect = sum(X == I);
        recover = sum(X == R);
        vaccine = sum(X == V);
        avers(t+1,:) = avers(t+1,:) + [suscept infect recover vaccine];
    end
end

% Average results
avers = avers/N;
newinfected = newinfected/N;
newvaccinated = newvaccinated/N;

% Plot results
close all;
plot(0:15,avers);
figure;
plot(0:15,newinfected);
figure;
plot(0:15,newvaccinated);

%%




%% 4 The H1N1 pandemic in Sweden 2009
S = 0;
I = 1;
R = 2;
V = 3;

Vacc = [5 9 16 24 32 40 47 54 59 60 60 60 60 60 60 60]'; % New Vacc-vector
I0 = [1 3 5 9 17 32 32 17 5 2 1 0 0 0 0 0]';    % Actual number of nodes who were infected each week.
N = 10; % Number of iterations during search.
n = 934;    % Number of nodes.
X0 = zeros(n,1);
setup = randi(n,round(n*Vacc(1)/100)+I0(1),1); % Find indices for intial state.
X0(setup(1)) = I;   % Choose patient zero.
X0(setup(2:end)) = V;   % Choose initially vaccinated individuals.
%%
kinit = 10;
b0 = 0.3;
r0 = 0.7;
param = [kinit b0 r0];  % Initial parameters for the gradient search.
lastparam = param;
% Initial parameter increments.
dk = 1;
db = 0.1;
dr = 0.1;
% Initial parameter spaces
kspace = param(1)-dk:dk:param(1)+dk;
bspace = param(2)-db:db:param(2)+db;
rspace = param(3)-dr:dr:param(3)+dr;
%%
precision = 7; % Choose how "precise" one wants to be.
RMSE = 10^4;
for loop = 1:precision
    RMSE = RMSE + 100;    % Increase RMSE so that we don't stop early on a fluke.
    notThere = 1;
    while(notThere) % We haven't found the ideal parameters for a given set of parameter increments.
        for i1 = 1:length(kspace)   % Creating graphs is a time-consuming business, minimize the number of times it's done.
            k = kspace(i1); % Create a graph for each k in kspace.
            
            c0 = k/2;
            k0 = round(k+1);
            W = ones(k0);
            verts = k0;
       
            for t = 1:n-k0
                if(mod(c0,1) ~= 0)
                    u = rand;
                    if(u < mod(c0,1))
                        c = ceil(c0);
                    else
                        c = floor(c0);
                    end
                else
                    c = c0;
                end
                
                W = [W zeros(size(W,2),1);
                     zeros(1,size(W,2)+1)];
                w = sum(W,2);
                P = w/sum(w);
                for j = 1:c
                    neighbor = randsample(1:verts+1, 1, true, full(P));
                    P(neighbor) = 0;
                    W(verts+1,neighbor) = 1;
                    W(neighbor,verts+1) = 1;
                end
                verts = verts + 1;
            end
            
            % For each combination of beta and rho, simulate the pandemic using the present random graph.
            for i2 = 1:length(bspace)
                for i3 = 1:length(rspace)
                    beta = bspace(i2); 
                    rho = rspace(i3);
        
                    G = sparse(W);
                    newinfected = zeros(16,1);
                    newinfected(1) = N;
                    for i = 1:N
                        X = X0;
                        for t = 1:15
                            infect = 0;
                            if(Vacc(t+1)-Vacc(t) > 0)
                                tovacc = round(Vacc(t+1)*n/100) - sum(X == V);
                                notvaccinated = find(X ~= V);
                                bigness = length(notvaccinated);
                                inds = notvaccinated(randi(bigness,tovacc,1));
                                X(inds) = V;
                            else
                                tovacc = 0;
                            end
        
                            Xtmp = (X == I);
                            m = G*Xtmp;
                            for q = 1:n
                                u = rand;
                                switch X(q)
                                    case 0
                                        if(u < (1-(1-beta)^(m(q))))
                                            X(q) = 1;
                                            infect = infect + 1;
                                        end
                                    case 1
                                        X(q) = (u < rho) + 1;
                                end
                            end
                                
                            newinfected(t+1) = newinfected(t+1) + infect;
                        end
                    end
                    newinfected = newinfected/N;
                    test = sqrt(sum((newinfected-I0).^2)/15);   % Calculate RMSE for present parameter set.
                    
                    % If it's the best set of parameters yet, store them
                    % and the resulting infection rate.
                    if(test < RMSE) 
                        param = [k beta rho];
                        It = newinfected;
                        RMSE = test;
                    end    
                end 
            end
        end
       
        if(sum(abs(param-lastparam)) == 0)
            notThere = 0;
        end
    
        kspace = param(1)-dk:dk:param(1)+dk;
        % Make sure that beta and rho stays in the interval [0,1] as they
        % are probabilities.
        if(param(2) >= 1)
            param(2) = 1;
            bspace = [param(2)-db,param(2)];
        elseif(param(2) + db > 1)
            bspace = [param(2)-db,param(2)];
        elseif(param(2) <= 0)
            param(2) = 0;
            bspace = [param(2),param(2)+db];
        elseif (param(2) - db < 0)
            bspace = [param(2),param(2)+db];
        else
            bspace = param(2)-db:db:param(2)+db;
        end
    
        if(param(3) >= 1)
            param(3) = 1;
            rspace = [param(3)-dr,param(3)];
        elseif(param(3) + dr > 1)
            rspace = [param(3)-dr,param(3)];
        elseif(param(3) <= 0)
            param(3) = 0;
            rspace = [param(3),param(3)+dr];
        elseif(param(3) - dr < 0)
            rspace = [param(3),param(3)+dr];
        else
            rspace = param(3)-dr:dr:param(3)+dr;
        end
        lastparam = param;
    end
    dk = dk/2;  % Calculate new parameter increments.
    db = db/2;
    dr = dr/2;
    switch loop     % For some iterations, store infection rates and current parameters.
        case 1
            I1 = It;
            param1 = param;
        case 3
            I3 = It;
            param3 = param;
        case 5
            I5 = It;
            param5 = param;
        case 7
            I7 = It;
            param7 = param;
    end
end

finalparam = param;
%% Plot best approximation result.
figure;
plot(0:15,It);
hold on;
plot(0:15,I0);

%% Final simulation
%simparam = [8.4688,0.2047,0.6063];  % Optimal parameters used for simulation.
simparam = [8.0059,0.1006,0.4002];
n = 934;
k = simparam(1);
c0 = k/2;
k0 = round(k+1);
W = ones(k0);
verts = k0;
 
for t = 1:n-k0
    if(mod(c0,1) ~= 0)
        u = rand;
        if(u < mod(c0,1))
            c = ceil(c0);
        else
            c = floor(c0);
        end
    else
        c = c0;
    end
       
    W = [W zeros(size(W,2),1);
        zeros(1,size(W,2)+1)];
    w = sum(W,2);
    P = w/sum(w);
    for j = 1:c
        neighbor = randsample(1:verts+1, 1, true, full(P));
        P(neighbor) = 0;
        W(verts+1,neighbor) = 1;
        W(neighbor,verts+1) = 1;
    end
    verts = verts + 1;
end
G = sparse(W);

%%
beta = simparam(2);
rho = simparam(3);       

N = 100;
setup = round(Vacc(1)*n/100)+I0(1);
inds = randi(n,setup,1);
X0 = zeros(n,1);
X0(inds(1)) = I;
X0(inds(2:end)) = V;

avers = zeros(16,4);
avers(1,:) = N*[sum(X0==S) sum(X0==I) sum(X0==R) sum(X0==V)];
newinfected = zeros(16,1);                    
newinfected(1) = N*sum(X==I);
newvaccinated = zeros(16,1);
for i = 1:N        
    X = X0;
    for t = 1:15
        infect = 0;
        if(Vacc(t+1)-Vacc(t) > 0)
            tovacc = round(Vacc(t+1)*n/100) - sum(X == V);
            notvaccinated = find(X ~= V);
            bigness = length(notvaccinated);
            inds = notvaccinated(randi(bigness,tovacc,1));
            X(inds) = V;    
        else
            tovacc = 0;    
        end
        
        Xtmp = (X == I);
        m = G*Xtmp;      
        for q = 1:n
            u = rand;
            switch X(q)
                case 0
                    if(u < (1-(1-beta)^(m(q))))
                        X(q) = 1;
                        infect = infect + 1;
                    end                    
                case 1                
                    X(q) = (u < rho) + 1;     
            end            
        end
        newinfected(t+1) = newinfected(t+1) + infect;
        newvaccinated(t+1) = newvaccinated(t+1) + tovacc;
        suscept = sum(X == S);
        infect = sum(X == I);
        recover = sum(X == R);
        vaccine = sum(X == V);
        avers(t+1,:) = avers(t+1,:) + [suscept infect recover vaccine];                   
    end    
end

avers = avers/N;
newinfected = newinfected/N;
newvaccinated = newvaccinated/N;

close all;
plot(0:15,avers);
title('Average SIRV distributions for best parameter estimate')
xlabel('Weeks');
ylabel('Number of individuals in subset');

figure;
plot(0:15,newinfected);
hold on;
plot(0:15,I0);
title('Average number of newly infected individuals vs real number')
xlabel('Weeks');
ylabel('Number of newly infected');
hold off;
