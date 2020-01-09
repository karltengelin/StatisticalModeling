load -ascii traffic.mat
load -ascii capacities.mat
load -ascii flow.mat
load -ascii traveltime.mat
load -ascii wifi.mat
%% Task 1a
weights = traveltime;

%we want to use graphshortestpath so therefore we need to create vector
%representations of where each link originate (which node) and where it
%points
for i = 1:size(traffic,2)          
    for j = 1:size(traffic,1)
        if traffic(j,i) == 1
            outlink(i) = j;
        elseif traffic(j,i) == -1
            inlink(i) = j;
        end
    end
end

%outlink shows outgoing links 
%inlink shows ingoing links

DG = sparse(outlink,inlink,weights,17,17);
[DIST,PATH,PRED] = graphshortestpath(DG,1,17)

%% Task 1b

%here we just use reformations of the data to siut the commando graphmaxflow

DG = sparse(outlink,inlink,capacities,17,17);
[M,F,K] = graphmaxflow(DG,1,17);
M

%% Task 1c
%here we are interested in the total flow in each node

ext_inflow = zeros(1,17);   %initializing external inflow i.e inflow into each node
ext_outflow = zeros(1,17);  %initializing external outflow i.e outflow out from each node

for i = 1:17
    X = find(i == inlink);
    Y = find(i == outlink);
    ext_inflow(i) = sum(flow(X));
    ext_outflow(i) = sum(flow(Y));
end

tot_flow = ext_outflow - ext_inflow;    %the total in-our outflow in each node

%% Task 1d - Compute social optimum

%first we identify our variables
B = traffic;
C = capacities;
nbroflinks = size(B,2);
l = traveltime;
lambda = zeros(length(ext_inflow),1);
mu = zeros(length(ext_outflow),1);

lambda(1) = tot_flow(1);
mu(end) = tot_flow(1);

%here we use the CVX environment to calculate optimal f_star

cvx_begin
variable f(nbroflinks)
minimize sum(l.*C.*inv_pos(1-f.*inv_pos(C)) - l.*C)
subject to
B*f == lambda - mu;
0 <= f <= C;
cvx_end
f_star = f;
f_star;
%% 1e - Compute Wardrop equilibrium
%again we use CVX to minimize the expression described in the report

cvx_begin
variable f(nbroflinks)
minimize sum(-C.*l.*log(1-f.*inv_pos(C)))
subject to
B*f == lambda - mu; 
0 <= f <= C;
cvx_end
f_zero = f;
f_zero;
%% 1f
%first we create our omega
omega = f_star.*l.*inv_pos(C).*(inv_pos(1-f_star.*inv_pos(C))).^2;      %omega = f_star.*(l.*C.*inv_pos((C - f).^2));

%Then we minimize using CVX

cvx_begin
variable f(nbroflinks)
minimize sum(-l.*C.*log(1-f.*inv_pos(C)) + f.*omega)  %sum(l.*C.*inv_pos(1-f.*inv_pos(C)) - l.*C + f.*omega)
subject to 
B*f == lambda - mu;
0 <= f <= C;
cvx_end
f
%% 1g
%first we need to calculate our f_star by minimzing the following
%expression
cvx_begin
variable f(nbroflinks)
minimize sum(l.*C.*inv_pos(1-f.*inv_pos(C)) - f.*l - l.*C)  
subject to 
B*f == lambda - mu;
0 <= f <= C;
cvx_end
f_star = f

%we create an omega given our calculated f_star
omega = f_star.*l.*C.*(inv_pos((C - f_star).^2));

cvx_begin
variable f(nbroflinks)
minimize sum(-l.*C.*log(1-f.*inv_pos(C)) - f.*l + f.*omega)  
subject to 
B*f == lambda - mu;
0 <= f <= C;
cvx_end
fw = f

norm(f_star-fw)