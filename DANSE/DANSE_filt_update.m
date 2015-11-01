function [node] = DANSE_filt_update(node,sim_param,DANSE_param,node_update)
% DANSE_FILT_UPDATE - update filter coeffcients of a node in batch mode
% Syntax:  [node] = DANSE_FILT_UPDATE(node,,sim_param,DANSE_param,node_update)
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
% Oct. 2015; Last revision: 31-Oct-2015
%------------- BEGIN CODE --------------
%% pre-allocate correlation matrix
Ryy = zeros(node(node_update).sensors+ ...
    (DANSE_param.nb_nodes-1)*DANSE_param.desired_sources,...
    node(node_update).sensors+...
    (DANSE_param.nb_nodes-1)*DANSE_param.desired_sources);
Rnn = Ryy;

% find all nodes ~= to idx_node
idx = find(node_update ~= 1:DANSE_param.nb_nodes);
% broadcast signals of other nodes in freqeuncy domain
%% update correlation maticies
Y = cat(1,node(node_update).ds_frame,cat(1,node(idx).ZY));
N = cat(1,node(node_update).n_frame,cat(1,node(idx).ZN));
for ll=1:sim_param.fftL/2+1
    Ryy = squeeze(Y(:,ll,:))*squeeze(Y(:,ll,:))';
    Rnn = squeeze(N(:,ll,:))*squeeze(N(:,ll,:))';
    
    Ryy = Ryy /  sim_param.ds_idx;
    Rnn = Rnn / sim_param.n_idx;
    % change this if estimating Rxx
    Rxx = Ryy;
    W_temp = (Rnn+Rxx)\ Rxx(:,1:DANSE_param.desired_sources);
    %local filters are the top entires of W_temp
    node(node_update).loc_filt_coeff(:,:,ll) = W_temp(1:node(node_update).sensors,:);
    %G coefficeints applied to other node z-signals
    Gkq_temp = W_temp(node(node_update).sensors+1:end,:);
    for jj = 1:numel(idx);
        node(node_update).gkq(idx(jj)).coeff(:,:,ll) = ...
            Gkq_temp((jj-1)*DANSE_param.desired_sources+1:jj*DANSE_param.desired_sources,:);
    end
%% Update broadcast signals
    for idx_z = 1:DANSE_param.desired_sources
        % find broadcast signals for every frame (sum(conj(W).*y = W'y)
        node(node_update).ZY(idx_z,ll,:) = ...
            sum(bsxfun(@times,conj(node(node_update).loc_filt_coeff(:,idx_z,ll)),node(node_update).ds_frame(:,ll,:)));
        node(node_update).ZN(idx_z,ll,:) = ...
            sum(bsxfun(@times,conj(node(node_update).loc_filt_coeff(:,idx_z,ll)),node(node_update).n_frame(:,ll,:)));
    end 
end
%------------- END OF CODE --------------