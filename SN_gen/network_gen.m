function [node,source,noise,white_noise_var] = ...
    network_gen(DANSE_param,sim_param,source_sig,noise_sig)
% NETWORK_GEN - generates random source, noise and node positions
% Syntax:  [node,source,noise,white_noise_var] = network_gen(DANSE_param,sim_param,source_sig,noise_sig)
% Inputs:   source_sig      -   desired source signals
%           noise_sig       -   additive noise (can be multichannel)
%           DANSE_param     -   structure with DANSE parameters
%           sim_param       -   structure with simulation parameters
%   
% Outputs:
%   node            -   contains the generated nodes in a structure
%   source          -   raw source signals
%   noise           -   raw noise signals
%   white_noise_var -   uncorrelated additive noise
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% May 2015; Last revision: 22-May-2015
%------------- BEGIN CODE --------------
% initalize node structure
node(sim_param.nb_nodes) = struct;
% generate sources at random locations 
for ii = 1:DANSE_param.desired_sources
    source(ii).nb = ii;                                 % source number
    source(ii).pos  = [sim_param.L(1)*rand(1) sim_param.L(2)*rand(1) 2];
    source(ii).signal = source_sig(:,ii); % source signal   
end

% generate noise sources at random locations
for ii = 1:DANSE_param.noise_sources
    noise(ii).nb = ii;                                  % noise source number
    noise(ii).pos  = [sim_param.L(1)*rand(1) sim_param.L(2)*rand(1) 2];
    noise(ii).signal =  noise_sig(:,ii);  % noise signal
end

white_noise_var = mean(var(cat(2,source.signal)))/2;    % representative of sensor noise (uncorrelated between nodes)

% generate nodes at random locations
for ii = 1:sim_param.nb_nodes
    node(ii).nb = ii;                                   % node number
    node(ii).sensors = DANSE_param.sensors;                               % number of sensors on node
    node(ii).ss_clean = zeros(sim_param.length_of_signal,node(ii).sensors); % pre-allocate clean source signals
    node(ii).ss_noise = zeros(sim_param.length_of_signal,node(ii).sensors); % pre-allocate additive noise signals
    node(ii).pos = [sim_param.L(1)*rand(1) sim_param.L(2)*rand(1) 2]; % node position
    
    % this is the sensor array for the node (can customize sensor positions) 
    r_sen = 0.1;                                    % sensor distance from center of node (equispaced) now 10 centimeters
    circ_pos = linspace(0,2*pi,node(ii).sensors+1); % position of sensors around node

    for jj = 1:node(ii).sensors
        node(ii).sensor(jj).pos(1) = node(ii).pos(1)+r_sen*cos(circ_pos(jj));   % x-axis sensor position on node
        node(ii).sensor(jj).pos(2) = node(ii).pos(2)+r_sen*sin(circ_pos(jj));   % y-axis sensor position on node
        node(ii).sensor(jj).pos(3) = 2;   % z-axis sensor position on node
        for kk = 1: size(source,2)
            d = norm(source(kk).pos - node(ii).sensor(jj).pos);                % calculate attenuation factor based on Euclidean distance
            node(ii).steering(jj,kk) = d;                                      % steering matrix (can be used to check G-coefficients in DANSE algorithm)
            node(ii).ss_clean(:,jj) = node(ii).ss_clean(:,jj)+d*source(kk).signal;
        end
        for kk = 1:size(noise,2)
            d = norm(noise(kk).pos - node(ii).sensor(jj).pos);
            node(ii).ss_noise(:,jj) = node(ii).ss_noise(:,jj)+d*noise(kk).signal+sqrt(white_noise_var)*randn(sim_param.length_of_signal,1);
        end
    end
    
    % generate local filter coeff and Gkqs for all other nodes
    node(ii).loc_filt_coeff = repmat(eye(node(ii).sensors,DANSE_param.desired_sources),1,1,sim_param.fftL/2+1);

    idx = find(ii ~= 1:DANSE_param.nb_nodes);
    for jj = idx
        node(ii).gkq(jj).coeff = zeros(DANSE_param.desired_sources,DANSE_param.desired_sources,sim_param.fftL/2+1);
    end 
end
%------------- END OF CODE --------------

