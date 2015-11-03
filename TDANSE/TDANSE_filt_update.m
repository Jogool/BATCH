function [node] = TDANSE_filt_update(node,sim_param,DANSE_param,node_update)
% TDANSE_FILT_UPDATE - update filter coeffcients of a node in batch mode 
% Syntax:  [node] = TDANSE_FILT_UPDATE(node,,sim_param,DANSE_param,node_update)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%           node_update     -   node to update
%                                                         
% Outputs:  node            -   return structure containing updated filter                                    coeffecieints and transmitted signals    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 03-Nov-2015
%------------- BEGIN CODE --------------
nb_ds = DANSE_param.desired_sources;
idx_nodes = node(node_update).tree_conn;
%% pre-allocate correlation matrix
Ryy = zeros(node(node_update).sensors+ ...
    numel(idx_nodes)*DANSE_param.desired_sources,...
    node(node_update).sensors+...
    numel(idx_nodes)*DANSE_param.desired_sources);
Rnn = Ryy;

% find all nodes connected to node (note since root node is updating, there
% is never any DF signals
%% update correlation maticies
Y = cat(1,node(node_update).ds_frame,cat(1,node(idx_nodes).ff_ZY));
N = cat(1,node(node_update).n_frame,cat(1,node(idx_nodes).ff_ZN));
for ll=1:sim_param.fftL/2+1
    Ryy = squeeze(Y(:,ll,:))*squeeze(Y(:,ll,:))';
    Rnn = squeeze(N(:,ll,:))*squeeze(N(:,ll,:))';
    
    Ryy = Ryy /  sim_param.ds_idx;
    Rnn = Rnn / sim_param.n_idx;
    % change this if estimating Rxx
    Rxx = Ryy;
    W_temp = (Rnn+Rxx)\ Rxx(:,1:nb_ds);
    %local filters are the top entires of W_temp
    node(node_update).loc_filt_coeff(:,:,ll) = W_temp(1:node(node_update).sensors,:);
    %G coefficeints applied to other node z-signals
    Gkq_temp = W_temp(node(node_update).sensors+1:end,:);
    for idx_node = 1:numel(idx_nodes);
        node(node_update).gkq(idx_nodes(idx_node)).coeff(:,:,ll) = ...
            Gkq_temp((idx_node-1)*nb_ds+1:idx_node*nb_ds,:);
    end
%% Update broadcast signals
    for idx_z = 1:nb_ds
        % find broadcast signals for every frame (sum(conj(W).*y = W'y)
        node(node_update).loc_ZY(idx_z,ll,:) = ...
            sum(bsxfun(@times,conj(node(node_update).loc_filt_coeff(:,idx_z,ll)),node(node_update).ds_frame(:,ll,:)));
        node(node_update).loc_ZN(idx_z,ll,:) = ...
            sum(bsxfun(@times,conj(node(node_update).loc_filt_coeff(:,idx_z,ll)),node(node_update).n_frame(:,ll,:)));
    end 
end
%------------- END OF CODE --------------