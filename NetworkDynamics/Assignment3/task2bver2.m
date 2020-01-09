load -ascii traffic.mat
load -ascii capacities.mat
load -ascii flow.mat
load -ascii traveltime.mat
load -ascii wifi.mat
load -ascii coords.mat
%%

iter = 1000;
W = wifi;

%set our inital states as random
states = zeros(100,iter);
states(:,1) = randi(8,100,1);

%initialize vectors and matrices
test = zeros(100,100);
sum1 = zeros(1,100);
sum2 = zeros(1,100);
sumsum = zeros(1,8);


colorcode = ['r' 'g' 'b' 'y' 'm' 'c' 'w' 'k'];


%% initial figure

    figure(1)
    set(gcf,'color','white')
    gplot(W,coords,'-k');
    hold on
        for iii = 1:100
            scatter(coords(iii,1),coords(iii,2),75,'markeredgecolor','k','markerfacecolor', colorcode(states(iii,1)));
        end
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    

%%


for i = 2:iter
    random_node = randi(100);
    states(:,i) = states(:,i-1);
    prob = zeros(length(colorcode),1);
    eta = -i/200;
    
    
    for color = 1:length(colorcode)
        for j = 1:100
            sum1(j) = W(random_node, j)*cost(color,states(j,i));
        end
        for jjj = 1:length(colorcode)
            for jj = 1:100
                sum2(jj) = W(random_node, jj)*cost(jjj,states(jj,i));
            end
            sumsum(jjj) = exp(eta.*sum(sum2));
        end
        prob(color) = (exp(eta.*sum(sum1)))/sum(sumsum);
    end
    
    F = cumsum(prob);
    u = rand;
    
    new_color = find(u<F,1);
    
    states(random_node,i) = new_color;
      
%         for iii = 1:size(states,1)
%             scatter(coords(iii,1),coords(iii,2),75,'markeredgecolor','k','markerfacecolor', colorcode(states(iii,i)));
%         end
        %pause(0.001);
    
    for k = 1:100
        for kk = 1:100
            test(k,kk)=W(k,kk).*cost(states(k,i),states(kk,i));
        end
    end
    
    U(i-1) = (1/2).*sum(sum(test));
    if U(i-1) == 0
        break
    end
    
end

for iii = 1:size(states,1)
    scatter(coords(iii,1),coords(iii,2),75,'markeredgecolor','k','markerfacecolor', colorcode(states(iii,end)));
end

figure
plot(U)

states(:,i)'



function answer = cost(a,X)
        if a == X
            answer = 2;
        elseif abs(X-a) == 1
            answer = 1;
        else
            answer = 0;
        end
end