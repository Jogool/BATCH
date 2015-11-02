function node = TIDANSE_init(node,sim_param,DANSE_param)
%TIDANSE_init - initilization of nodes for TI-DANSE algorithm 
% Syntax:  [node] = TIDANSE_init(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%                                                         
% Outputs:  node            -   return structure containing initial cost 
%                               and transmitted signals                        
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 02-Nov-2015
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
                sum(bsxfun(@times,conj(node(idx_node).P(:,idx_z,ll)),node(idx_node).ds_frame(:,ll,:)));
            node(idx_node).loc_ZN(idx_z,ll,:) = ...
                sum(bsxfun(@times,conj(node(idx_node).P(:,idx_z,ll)),node(idx_node).n_frame(:,ll,:)));
        end
    end
    % resest G-coefficents for TI-DANSE
    node(idx_node).gkq = [];
    node(idx_node).gkq(1).coeff = zeros(DANSE_param.desired_sources,DANSE_param.desired_sources,sim_param.fftL/2+1);
end
% calculate cost
node = TIDANSE_cost(node,sim_param,DANSE_param);
%------------- END OF CODE --------------