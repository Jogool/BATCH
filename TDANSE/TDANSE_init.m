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
nb_ds = DANSE_param.desired_sources;
%% generate initial local broadcast signals
for idx_node = 1:DANSE_param.nb_nodes;    
    % initalize broadcast signals
    node(idx_node).loc_ZY = zeros(nb_ds,sim_param.fftL/2+1,sim_param.ds_idx);
    node(idx_node).loc_ZN = zeros(nb_ds,sim_param.fftL/2+1,sim_param.n_idx);
    for idx_z = 1:nb_ds
        for ll=1:sim_param.fftL/2+1
            % find broadcast signals for every frame (sum(conj(W).*y = W'y)
            node(idx_node).loc_ZY(idx_z,ll,:) = ...
                sum(bsxfun(@times,conj(node(idx_node).loc_filt_coeff(:,idx_z,ll)),node(idx_node).ds_frame(:,ll,:)));
            node(idx_node).loc_ZN(idx_z,ll,:) = ...
                sum(bsxfun(@times,conj(node(idx_node).loc_filt_coeff(:,idx_z,ll)),node(idx_node).n_frame(:,ll,:)));
        end
    end
end
node = TDANSE_cost(node,sim_param,DANSE_param);
%------------- END OF CODE --------------
