function [] = plot_WSN(node,A,source,noise)
% PLOT_WSN - plot generated network
% Syntax:  [] = PLOT_WSN(node,A,source,noise)
% Inputs:   node            -   contains the generated nodes in a structure
%           A               -   adjacency matrix of tree
%           source          -   desired sources in a structure
%           noise           -   noise sources in a structure
%   
% Outputs:  none
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
clf
axis([0 5 0 5]);
hold on
A = tril(A);
% plot desired sources
for ii = 1:size(source,2)
    plot(source(ii).pos(1),source(ii).pos(2),'s','MarkerFaceColor','red', ...
        'MarkerEdgeColor','red')
end

% plot noise sources
for ii = 1:size(noise,2)
    plot(noise(ii).pos(1),noise(ii).pos(2),'d','MarkerFaceColor','blue', ...
        'MarkerEdgeColor','blue')
end

% plot nodes
for ii = 1:size(node,2)
    plot(node(ii).pos(1),node(ii).pos(2),'o','MarkerFaceColor','black', ...
        'MarkerEdgeColor','black')
    text(node(ii).pos(1)+0.1,node(ii).pos(2)+0.1,num2str(ii))
end

% plot ad-hoc connections
for ii = 1:size(node,2)
    nb_conn = find(A(ii,:));
    for jj = nb_conn
        h = line([node(ii).pos(1) node(jj).pos(1)],[node(ii).pos(2) node(jj).pos(2)]);
        set(h,'Color','red')
        set(h,'LineWidth',1)
        set(h,'LineStyle','--')
    end

% plot tree connections
    nb_conn = node(ii).tree_conn';
    for jj = nb_conn
        h = line([node(ii).pos(1) node(jj).pos(1)],[node(ii).pos(2) node(jj).pos(2)]);
        set(h,'Color','black')
        set(h,'LineWidth',2)
    end
    
end
drawnow
%------------- END OF CODE --------------
