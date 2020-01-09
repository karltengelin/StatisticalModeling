clear all
close all
load('IOdownload.mat')
%% 1.a The in-degree and out-degree centrality
% First we summarize the number of out and in degrees for each sector for each country
in_swe = sum(io.swe2000,1);
out_swe = sum(io.swe2000,2)';
in_ind = sum(io.idn2000,1);
out_ind = sum(io.idn2000,2)'; 

% Secondly we want to check which three nodes that have the highest out
% degree and in degree
[out_max_swe,out_max_pos_swe] = maxk(out_swe,3);
[in_max_swe,in_max_pos_swe] = maxk(in_swe,3);
[out_max_ind,out_max_pos_ind] = maxk(out_ind,3);
[in_max_ind,in_max_pos_ind] = maxk(in_ind,3);

%Then we check which name that corresponds to the selected nodes and store in
%a name-vector
for i = 1:3
name_in_swe(i) = name(in_max_pos_swe(i));
name_out_swe(i) = name(out_max_pos_swe(i));
name_in_ind(i) = name(in_max_pos_ind(i));
name_out_ind(i) = name(out_max_pos_ind(i));
end

%lastly we display the results
disp('in centrality swe:')
disp(name_in_swe)
disp('out centrality swe:')
disp(name_out_swe)
disp('in centrality ind:')
disp(name_in_ind)
disp('out centrality ind:')
disp(name_out_ind)
%% 1.b The eigenvector centrality on the largest connected component;

% First we want to make sure that we work with the largest connected
% component, and by looking in the W matrix for sweden we can see that there are a
% couple of nodes (10 for example) that isn't connected to any other.
%so the first step is to remove these nodes 
swe2000_new = io.swe2000;
swe2000_new(:,all(~swe2000_new,1))=[];
swe2000_new(all(~swe2000_new,2),:)=[];

%By looking in the W matrix for indonesia we can see that the graph is
%strongly connected. This can be done by looking at the three first nodes
%for a start, here we see that they all point to each other and if we look
%further down in the matrix the only node that does not point to any of
%them is node 21. However node 21 is connected to several nodes (node 8 for
%example) which in turn is connected to node 1, 2 and 3, thus the graph is
%stronly connected.
idn2000_new = io.idn2000;

% Since we now have shuffled around the nodes in the swedish W we need to change the name
% list so that the correct corresponds to the correct name:
name_new = name;
name_new(sum(io.swe2000,2)==0) = [];

% now we are interested in the eigenvalue of each W, for this we use the
% command eig:
[V1,D1] = eig(swe2000_new');
[V2,D2] = eig(idn2000_new');
% we now want to check the leading eigenvalues (which is real valued), and since our D-matrix is
% a diagonal matrix the command sum will give us the eigenvalue for each
% row. And since we are interested in the leading eigenvalue we check which
% row (i.e which index) that is the biggest 
[value1 index1] = max(sum(abs(D1),2));
[value2 index2] = max(sum(abs(D2),2));
%store the eigenvectors assosiated with the calculated eigenvalues in a
%vector
eigenvect1 = abs(V1(:,index1));
eigenvect2 = abs(V2(:,index2));
%taking out the three biggest values in the vectors according to the
%eigenvector centrality algorithm
[vectval1,vectpos1] = maxk(eigenvect1,3);
[vectval2,vectpos2] = maxk(eigenvect2,3);

%Connect the indexes with the correct names:
for i = 1:3
name_swe(i) = name_new(vectpos1(i));
name_ind(i) = name(vectpos2(i));
end

%displaying the result:
disp('eigenvectcentrality swe:')
disp(name_swe)
disp('eigenvectcentrality indoniesia:')
disp(name_ind)
%% 1.c The Katz centrality, with ? = 0.15 and with two different values on ?, 
%namely ? = 1 and when ?i = 1 for the ?Wholesale & retail trade; repairs? sector and 
%zero for all other sectors.

beta = 0.15;
lambda_swe = value1; %here we take the eigenvalues from the previous task
lambda_idn = value2;
I = eye(size(io.swe2000));
mu_1 = ones(1,length(io.swe2000))'; % first mu is a vector of ones
mu_2 = zeros(1,length(io.swe2000))'; 
mu_2(31) = 1; % second mu is a vector with zeros everywhere except on place 31 (?Wholesale & retail trade; repairs?)

%now we create the z vectors which contains our centrality-information, the
%equations follow equation 2.3 in lecture notes
z_b_1_swe = (I-lambda_swe^-1*(1-beta)*io.swe2000')\(mu_1*beta); 
z_b_2_swe = (I-lambda_swe^-1*(1-beta)*io.swe2000')\(mu_2*beta);

z_b_1_idn = (I-lambda_idn^-1*(1-beta)*io.idn2000')\(mu_1*beta);
z_b_2_idn = (I-lambda_idn^-1*(1-beta)*io.idn2000')\(mu_2*beta);

%We once again extract the three largest values
[max_mu_1_val_swe max_mu_1_pos_swe] = maxk(z_b_1_swe,3);
[max_mu_2_val_swe max_mu_2_pos_swe] = maxk(z_b_2_swe,3);

[max_mu_1_val_idn max_mu_1_pos_idn] = maxk(z_b_1_idn,3);
[max_mu_2_val_idn max_mu_2_pos_idn] = maxk(z_b_2_idn,3);

%connecting the index to a name in the name vector
for i = 1:3
name_swe_mu1(i) = name(max_mu_1_pos_swe(i));
name_swe_mu2(i) = name(max_mu_2_pos_swe(i));
name_idn_mu1(i) = name(max_mu_1_pos_idn(i));
name_idn_mu2(i) = name(max_mu_2_pos_idn(i));
end

%displaying the result
disp('katz swe mu_1:')
disp(name_swe_mu1)
disp('katz swe mu_2:')
disp(name_swe_mu2)
disp('katz ind mu_1:')
disp(name_idn_mu1)
disp('katz ind mu_2:')
disp(name_idn_mu2)