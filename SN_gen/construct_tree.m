function [node,A,uo] = construct_tree(node)
% CONSTRUCT_TREE - construct a MST based on current node positions 
% Syntax:  [node,A,uo] = construct_tree(node)
% Inputs:   node            -   contains the generated nodes in a structure
%                                                         
% Outputs:
%           node            -   returns node with updated tree connections 
%                                   in structure format
%           A               -   adjacency matrix of tree
%           uo              -   update order for TDANSE
%
% Other m-files required: path_find,graphminspantree,sparse
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 04-Dec-2014
%------------- BEGIN CODE --------------
nb_nodes = size(node,2);
r = 0;                              % initial broadcast radius of each node
alg_conn = 0;                       % algebratic connectivity
while lt(alg_conn,1e-10)

    A = zeros(nb_nodes);            % adjacency matrix based on current connectivity      
    D = zeros(nb_nodes);            % matrix that holds the distances between nodes
    for ii = 1:nb_nodes
        index = find(ii ~= 1:nb_nodes);
        for jj = 1:nb_nodes-1
            D(ii,index(jj)) = norm(node(ii).pos-node(index(jj)).pos);
            if lt(D(ii,index(jj)),r)    % if the distance is less than the sensing radius, assume nodes can communicate
                A(ii,index(jj)) = 1;
            end
        end
    end

    Deg_mat = diag(sum(A,2));   % diagnomal matrix is equal to connections
    L = Deg_mat-A;              % Laplacian matrix 
    lambda = eig(L);            % eigenvalues of Laplacian matrix
    [~,I] = sort(lambda);
    alg_conn = lambda(I(2));    % algebratic connectivity is equal to the second largest eigenvalue of Laplacian matrix
    r = r + 0.01;               % increase sensing radius 
end

% reduce node communication range to minimal distance for connectivity
max_range = zeros(nb_nodes,1);
for ii = 1:nb_nodes
    index = find(A(ii,:));
    for jj = 1:length(index)
        max_range(ii) = max(max_range(ii),  D(ii,index(jj)));
    end
end

%find all posssible connections (ad-hoc)
for ii = 1:nb_nodes
    node(ii).conn = find(A(ii,:));
end

% find miminum spanning tree based on current connected adjanceny matrix
% Euclidean weights
A_tril = tril(A);
euc_weight = D(find(A_tril));
[rows,cols] = find(A_tril);

UG = sparse(rows,cols,euc_weight,nb_nodes,nb_nodes);

[ST,~] = graphminspantree(UG);

A_mst = full(ST)+full(ST)';    % adjanceny matrix of tree topology
A_mst(find(A_mst)) = 1;

%find all connections for tree
for ii = 1:nb_nodes
    node(ii).tree_conn = find(A_mst(:,ii));
end

[u1,v1]=eig(A_mst);
[~,lambdaind]=sort(abs(diag(v1)));
EigvCentrality=[[1:size(A_mst,1)]' abs(u1(:,lambdaind(end)))];

[~,root_node] = max(EigvCentrality(:,2));   % pick root node based on eigenvalue centrality

% find updating order based on root node and 
uo = root_node;
uo = path_find(uo,A_mst,root_node);
end
%------------- END OF CODE --------------

