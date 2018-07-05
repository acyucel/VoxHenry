function [nodeid_4_grnd,nodeid_4_injectcurr]=lse_assign_exc_grnd_nodes(nodeid_lft,nodeid_rght,nodeid_wlcond,num_ports,port_no)


nodeid_4_grnd=[]; nodeid_4_injectcurr=[];

nodeid_4_injectcurr=nodeid_lft{port_no}(:);

for kk=1:num_ports
    nodeid_4_grnd=[nodeid_4_grnd;nodeid_rght{kk}(:)];
end

for kk=1:num_ports
    if (kk ~= port_no)
        nodeid_4_grnd=[nodeid_4_grnd;nodeid_lft{kk}(:)];
    end
end

if (isempty(nodeid_wlcond) == 0)
    nodeid_4_grnd=[nodeid_4_grnd;nodeid_wlcond];
end