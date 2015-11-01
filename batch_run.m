function [] = batch_run(source_sig,noise_sig)
% BATCH_RUN - wrapper for batch mode simulations 
% Syntax:  [] = batch_run(source_sig,noise_sig)
% Inputs:   source_sig      -   desired signal
%           noise_sig       -   noise signal
%                                                         
% Outputs:  node            -   none
%
% Other m-files required: gen_param,network_gen,construct_tree,plot_WSN,
%                           centralized_batch,DANSE_init,DANSE_batch,
%                           TDANSE_init,TDANSE_batch
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% May. 2015; Last revision: 11-Oct-2015
%------------- BEGIN CODE --------------
% generate parameters for simulation environment and nodes
[sim_param, DANSE_param] = gen_param(source_sig,noise_sig);

% generate random network, construct tree and define updating order to
% T-DANSE algorithm
if sim_param.plot
    close all
    [node,source,noise,wnv] = ...
        network_gen(DANSE_param,sim_param,source_sig,noise_sig);
    [node,A,updateorder] = construct_tree(node);
    plot_WSN(node,A,source,noise)
else
    [node,~,~,~] = network_gen(DANSE_param);
    [node,~,updateorder] = construct_tree(node);
end

[node,sim_param] = frame_signal(node,sim_param);
%% Centralized
disp('Centralized')
node = centralized_batch(node,DANSE_param,sim_param);
cent_cost = [node.cent_cost];
sum_cent_cost = sum([node.cent_cost]);
% save orignal node
org_node = node;
%% DANSE
fprintf('\n')
disp('DANSE')
reverseStr = '';
% initialize nodes for DANSE algorithm
[node] = DANSE_init(node,sim_param,DANSE_param);
% inital cost
cost_DANSE = cat(1,node.cost);
% inital sum of costs
cost_sum_DANSE = sum(cost_DANSE);
% set node to start updating
node_update = updateorder(1);
ii = 1;
% set convergence flag
convergence = 0;
% check if we have met either condition
while ~or(convergence,ge(ii,DANSE_param.max_iter))
    [node] = DANSE_batch(node,sim_param,DANSE_param,node_update);
    cost_DANSE = [cost_DANSE cat(1,node.cost)];
    cost_sum_DANSE = [cost_sum_DANSE sum(cat(1,node.cost))];
    convergence =  norm(cat(1,node.cent_cost) - ...
        cellfun(@(x) x(end), {node.cost})') < DANSE_param.thresh;
    ii = ii + 1;  
    % round robin updating
    node_update=rem(node_update,sim_param.nb_nodes)+1;   
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n')
% reset node
node = org_node;
%% T-DANSE
fprintf('\n')
disp('T-DANSE')
reverseStr = '';
% initialize nodes for TDANSE algorithm
[node] = TDANSE_init(node,updateorder(1),sim_param,DANSE_param);
% inital cost
cost_TDANSE = cat(1,node.cost);
% inital sum of costs
cost_sum_TDANSE = sum(cost_TDANSE);
% set node to start updating
node_update = updateorder(1);
ii = 1;
% set convergence flag
convergence = 0;
% check if we have met either condition
while ~or(convergence,ge(ii,DANSE_param.max_iter))
    [node] = TDANSE_batch(node,sim_param,DANSE_param,node_update);
    cost_TDANSE = [cost_TDANSE cat(1,node.cost)];
    cost_sum_TDANSE = [cost_sum_TDANSE sum(cat(1,node.cost))];
    convergence =  norm(cat(1,node.cent_cost) - ...
        cellfun(@(x) x(end), {node.cost})') < DANSE_param.thresh;
    ii = ii + 1;  
    % follow path of tree for updating order
    node_update=updateorder(rem(ii,numel(updateorder))+1);
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n')
%% End matter
figure
hold on
loglog(cost_sum_DANSE)
loglog(cost_sum_TDANSE,'-xm')
axis tight
h =  refline(0,sum_cent_cost);
set(h,'LineStyle','--');

a = get(gca,'YLim');
set(gca,'YLim',[sum([node.cent_cost]) - sum([node.cent_cost])*.1 a(2)])
legend('DANSE', 'T-DANSE','Optimal');
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Iteration')
ylabel('Sum of LS cost for all nodes (dB)')
axis tight
%------------- END OF CODE --------------

