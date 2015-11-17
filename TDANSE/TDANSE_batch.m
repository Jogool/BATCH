function node = TDANSE_batch(node,sim_param,DANSE_param,node_update)
% TDANSE_BATCH - Wrapper to run a single iteration of the T-DANSE algorithm 
% Syntax:  node = TDANSE_batch(node,sim_param,DANSE_param,node_update)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%           node_update     -   node to update
%                                                         
% Outputs:  node            -   return structure containing cost, updated
%                                   filter and broadcast signals
%
% Other m-files required: TDANSE_ff, TDANSE_filt_update, TDANSE_cost
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 17-Nov2015
%------------- BEGIN CODE --------------
% fusion flow toward updating (root) node
node = TDANSE_ff(node,node_update,sim_param,DANSE_param);
% filter update
node = TDANSE_filt_update(node,sim_param,DANSE_param,node_update);
% cost calculation
node = TDANSE_cost(node,sim_param,DANSE_param);
%------------- END OF CODE --------------
