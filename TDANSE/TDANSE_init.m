function node = TDANSE_init(node,root,sim_param,DANSE_param)
% TDANSE_init - initilization of nodes for T-DANSE algorithm 
% Syntax:  [node] = TDANSE_init(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           root            -   root (sink) node of the tree
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%                                                         
% Outputs:  node            -   return structure containing initial cost 
%                               and transmitted signals                        
%
% Other m-files required: filterwithW
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 12-Oct-2015
%------------- BEGIN CODE --------------
% inital cost variable
[node.cost] = deal(0);
%% generate initial local broadcast signals
for idx_node = 1:DANSE_param.nb_nodes;    
    for frame_idx = 1:sim_param.ds_idx;
        for ll=1:sim_param.fftL/2+1
            node(idx_node).loc_ZY(:,ll,frame_idx) = ...
                node(idx_node).loc_filt_coeff(:,:,ll)'*node(idx_node).ds_frame(:,ll,frame_idx);
        end
    end
    for frame_idx = 1:sim_param.n_idx;
        for ll=1:sim_param.fftL/2+1
            node(idx_node).loc_ZN(:,ll,frame_idx) = ...
                node(idx_node).loc_filt_coeff(:,:,ll)'*node(idx_node).n_frame(:,ll,frame_idx);
        end
    end
    % local broadcast signals in time domain (only used for cost calculation)
    node(idx_node).loc_zy = ...
        filterwithW([node(idx_node).ss_clean],node(idx_node).loc_filt_coeff);
    node(idx_node).loc_zn = ...
        filterwithW([node(idx_node).ss_noise],node(idx_node).loc_filt_coeff);
    
end
node = TDANSE_ff(node,root,sim_param,DANSE_param);
node = TDANSE_df(node,root,sim_param);

for idx_node = 1:DANSE_param.nb_nodes;  
    idx = node(idx_node).ff_rec;
    % gather local signals and z-signals of other nodes
    y = [node(idx_node).ss_clean + node(idx_node).ss_noise];
    zy = [node(node(idx_node).ff_trans).df(idx_node).zy node(idx).ff_zy];
    zn = [node(node(idx_node).ff_trans).df(idx_node).zn node(idx).ff_zn];
    
    % gather G coefficiens
    Gkq_coeff = ...
        cat(1,[node(idx_node).gkq(node(idx_node).ff_trans).coeff],cat(1,node(idx_node).gkq(idx).coeff));
    
    W = cat(1, node(idx_node).loc_filt_coeff, Gkq_coeff);
    % estimated desired signal
    dest = filterwithW([y zy+zn],W);
    % calculate cost
    node(idx_node).cost = ...
        norm(node(idx_node).ss_clean(:,1:DANSE_param.desired_sources) - dest)^2;
end
%------------- END OF CODE --------------
