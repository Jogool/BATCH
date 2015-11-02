function [node] = TIDANSE_filt_update(node,sim_param,DANSE_param,node_update)
% TIDANSE_FILT_UPDATE - update filter coeffcients of a node in batch mode 
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
% Oct. 2015; Last revision: 02-Nov-2015
%------------- BEGIN CODE --------------
%% pre-allocate correlation matrix
Ryy = zeros(node(node_update).sensors+ DANSE_param.desired_sources,...
        node(node_update).sensors+DANSE_param.desired_sources);
Rnn = Ryy;

% index of all other non-updating nodes
idx = find(node_update ~= 1:DANSE_param.nb_nodes);

%% collect broadcast signals from other nodes and compute the sum
ZY = cat(1,node(idx).loc_ZY);
ZN = cat(1,node(idx).loc_ZN);
for idx_DANSE = 1:DANSE_param.desired_sources;
    ZY_sum(idx_DANSE,:,:) = sum(ZY(idx_DANSE:DANSE_param.desired_sources:end,:,:),1);
    ZN_sum(idx_DANSE,:,:) = sum(ZN(idx_DANSE:DANSE_param.desired_sources:end,:,:),1);
end
%% update correlation maticies
Y = cat(1,node(node_update).ds_frame,ZY_sum);
N = cat(1,node(node_update).n_frame,ZN_sum);

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
    node(node_update).gkq(1).coeff(:,:,ll) = W_temp(node(node_update).sensors+1:end,:);
    node(node_update).P(:,:,ll) = node(node_update).loc_filt_coeff(:,:,ll)/node(node_update).gkq(1).coeff(:,:,ll);
    for idx_z = 1:DANSE_param.desired_sources
        % find broadcast signals for every frame (sum(conj(W).*y = W'y)
        node(node_update).loc_ZY(idx_z,ll,:) = ...
            sum(bsxfun(@times,conj(node(node_update).P(:,idx_z,ll)),node(node_update).ds_frame(:,ll,:)));
        node(node_update).loc_ZN(idx_z,ll,:) = ...
            sum(bsxfun(@times,conj(node(node_update).P(:,idx_z,ll)),node(node_update).n_frame(:,ll,:)));
    end
end
% %% update filter 
% for ll=1:sim_param.fftL/2+1
%     % update local filter coefficients
%     W_temp = (Rnn(:,:,ll)+Rxx(:,:,ll))\...
%         Rxx(:,1:DANSE_param.desired_sources,ll);
%     % local filters are the top entires of W_temp
%     node(node_update).loc_filt_coeff(:,:,ll) = W_temp(1:node(node_update).sensors,:);
%     % G coefficients applied to sum of z-signals
%     node(node_update).gkq(1).coeff(:,:,ll) = W_temp(node(node_update).sensors+1:end,:);
%     % broadcast coefficients for TI-DANSE
%     node(node_update).P(:,:,ll) = node(node_update).loc_filt_coeff(:,:,ll)/node(node_update).gkq(1).coeff(:,:,ll);
%     
%     %% Update broadcast signals
%     for idx_z = 1:DANSE_param.desired_sources
%         % find broadcast signals for every frame (sum(conj(W).*y = W'y)
%         node(node_update).loc_ZY(idx_z,ll,:) = ...
%             sum(bsxfun(@times,conj(node(node_update).P(:,idx_z,ll)),node(node_update).ds_frame(:,ll,:)));
%         node(node_update).loc_ZN(idx_z,ll,:) = ...
%             sum(bsxfun(@times,conj(node(node_update).P(:,idx_z,ll)),node(node_update).n_frame(:,ll,:)));
%     end 
% end
% %% update broadcast signal in frequency domain
% for frame_idx = 1:sim_param.ds_idx;
%     for ll=1:sim_param.fftL/2+1
%         node(node_update).loc_ZY(:,ll,frame_idx)  = ...
%             node(node_update).P(:,:,ll)'*node(node_update).ds_frame(:,ll,frame_idx);
%     end
% end
% for frame_idx = 1:sim_param.n_idx;
%     for ll=1:sim_param.fftL/2+1
%         node(node_update).loc_ZN(:,ll,frame_idx)  = ...
%             node(node_update).P(:,:,ll)'*node(node_update).n_frame(:,ll,frame_idx);
%     end
% end
%------------- END OF CODE --------------