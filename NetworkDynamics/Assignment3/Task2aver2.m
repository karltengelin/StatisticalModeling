load -ascii traffic.mat
load -ascii capacities.mat
load -ascii flow.mat
load -ascii traveltime.mat
load -ascii wifi.mat
load -ascii coords.mat
%% 
% we initialize our variables
iter = 10000;
W = zeros(10);
v = ones(1,9);
D1 = diag(v,1);
D2 = diag(v,-1);
states = zeros(10,iter);
states(:,1) = ones(10,1);

%we create our toeplitz matrix
W = W+D1+D2;

%we specify our color set
colors = [1:2]; %1-red, 2-green

cordsX = (1:10)';
cordsY = zeros(10,1);
cords = [cordsX cordsY];

%% initial figure

    figure(1)
    set(gcf,'color','white')
    gplot(W,cords,'-k');
    hold on
        for iii = 1:10
            if states(iii,1) == 1
                scatter(cords(iii,1),cords(iii,2),200,'markeredgecolor','k','markerfacecolor', 'r');
            else
                scatter(cords(iii,1),cords(iii,2),200,'markeredgecolor','k','markerfacecolor', 'g');
            end

        end
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    

%%

%here we loop over pre-specified number of iterations
for i = 2:iter
    random_node = randi(10);        %create a radnom variable which selects a random node
    states(:,i) = states(:,i-1);    %reloading the states for each time step
    prob = zeros(length(colors),1); %the probability vector which shows the probability for a node to take on each color
    eta = -i/100;                   %our eta function
    
    %Here we create the probability that the node transitions to each given
    %color given the costfunction (here 'color == states(j,i)' which is a
    %one if they are the same otherwise zero, and the weight matrix (see the report for further details). 
    %we then store in our probability vector ''prob''
    for color = 1:length(colors)
        for j = 1:10
            sum1(j) = W(random_node, j)*(color == states(j,i));
        end
        for jjj = 1:length(colors)
            for jj = 1:10
                sum2(jj) = W(random_node, jj)*(jjj == states(jj,i));
            end
            sumsum(jjj) = exp(eta.*sum(sum2));
        end
        prob(color) = (exp(eta.*sum(sum1)))/sum(sumsum);
    end
    
    
    %we create a CDF out of the prob vector
    F = cumsum(prob);
    %we draw u from a uniform distribution
    u = rand;
    
    %the new color of the random node is given by:
    new_color = find(u<F,1);
    %we update the random node
    states(random_node,i) = new_color;
              
                for iii = 1:10              %we simulate the process
                    if states(iii,i) == 1
                        scatter(cords(iii,1),cords(iii,2),200,'markeredgecolor','k','markerfacecolor', 'r');
                    else
                        scatter(cords(iii,1),cords(iii,2),200,'markeredgecolor','k','markerfacecolor', 'g');
                    end

                end
                pause(0.2);
        %here we create the parts we need to create our potential function,
        %we store each value in a matrix in which we later summarize
      for k = 1:10
        for kk = 1:10
            test(k,kk)=W(k,kk).*cost(states(k,i),states(kk,i));
        end
      end
    
    %we create our potential function
    U(i-1) = (1/2).*sum(sum(test));
    if U(i-1) == 0      %if the potential function reaches zero we break
        break
    end
end
figure
plot(U)

states(:,i)'

%our cost function
function answer = cost(a,X)
        if a == X
            answer = 1;
        elseif a ~= X
            answer = 0;
        end
end