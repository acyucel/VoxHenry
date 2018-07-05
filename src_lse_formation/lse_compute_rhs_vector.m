function [rhs_vect] = lse_compute_rhs_vector(Ae,nodeid_4_injectcurr)

% Define right hand side
num_node = size(Ae,1);
num_curr = size(Ae,2);
rhs_vect=zeros(num_curr+num_node,1);

% sum the entries of excitation nodes and assign to rhs
tmp_mat=-Ae';
exc_mat=tmp_mat(:,nodeid_4_injectcurr);
rhs_vect(1:num_curr)=full(sum(exc_mat,2));
clear tmp_mat
