function [currs_port_yparams2] = lse_compute_Y_mat_column_alternative(num_ports,Ae,Ae_only_leaving,Ae_only_entering_bndry,x,nodeid_lft,currs_port_yparams)

num_node = size(Ae,1); num_curr = size(Ae,2);
x_inner_node=Ae_only_leaving*x(1:num_curr);
x_bndry_node=Ae_only_entering_bndry*x(1:num_curr);
x_node=-x_inner_node+x_bndry_node;

currs_port_yparams2=zeros(num_ports,1);

for kk=1:num_ports
    currs_port_yparams2(kk,1)=sum(x_node(nodeid_lft{kk}));
end

if (abs(currs_port_yparams2-currs_port_yparams(:,1)) > 1e-12)
    diff=abs(currs_port_yparams2-currs_port_yparams(:,1));
    %disp(['Check this difference :::',diff,'why? - nodes in nodeidlft should be excitation'])
end