Alpha = [0 2/5 1/5 0 0; 0 0 3/4 1/4 0; 1/2 0 0 1/2 0; 0 0 1/3 0 2/3; 0 1/3 0 1/3 0];
w = Alpha*ones(size(Alpha,2),1);
r = max(w);
D = diag(w);
P_targets = D\Alpha;

% number of iterations or expected running time
Tmax = 60; 


% Transition probability matrix

P_transprob = zeros(size(P_targets));
for i = 1:length(P_transprob)
    for j = 1:length(P_transprob)
        if i ~= j
            P_transprob(i,j) = Alpha(i,j)/max(w);
        end
    end
    P_transprob(i,i) = 1-sum(P_transprob(i,:));
end

%% --- Start of simulations --- %%

nbr_of_part = 100;
nbr_of_sim = 1000;
nodes = 1:length(Alpha);                
nbr_in_node = zeros(length(Alpha),Tmax);
start = 1;                                  %specifying starting node: 1-o, 2-a, 3-b, 4-c, 5-d                            %which is our end node? 1-o, 2-a, 3-b, 4-c, 5-d

for iii = 1:nbr_of_sim          
    for ii = 1:nbr_of_part
    %---- Simulate the particle moving around ----%

    x = zeros(5,Tmax);
    x(start,1) = 1;         
    n = start;         

        for i = 2:Tmax
            random = rand;
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
            x(n, i) = 1;
        end
        nbr_in_node = nbr_in_node + x;
    end
    
            % If we are interested in only plotting the number of particles in 
            % each node for a single run (i.e notan average), we only run the inner
            % loop and use this plotting:
% figure                                    
% hold on
%     for i = 1:5
%         plot(nbr_in_node(i,:))
%     end
% legend('node o', 'node a', 'node b', 'node c', 'node d','FontSize',12)
% ylabel('nbr of particles')
% xlabel('time step')
% hold off

end

Q_end = nbr_in_node(:,end)./nbr_of_sim;         % Gives the average number of particles in each node AT THE END OF 60 time steps
Q_full = nbr_in_node./nbr_of_sim;               % Gives the average number of particles in each node DURING 60 time steps
average_in_each_node_at_end = [nodes; Q_end']   % Presentation of results
                                                
                                                    
figure                                          % plots the average number of particles in each node during 60 time steps
hold on
    for i = 1:5
        plot(Q_full(i,:))
    end
legend('node o', 'node a', 'node b', 'node c', 'node d','FontSize',12)
ylabel('nbr of particles')
xlabel('time step')
hold off



