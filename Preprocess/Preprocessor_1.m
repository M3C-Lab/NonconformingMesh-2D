function [DataArrays, BCs] = Preprocessor_1(msh)
% Preprocessor for 1-domain

fprintf("\nPreprocess routine:\n")

% Degree of freedom at node
DataArrays.dof = 1;

% Number of nodes
DataArrays.nNode = msh.nbNod;
fprintf("  Total nodes: %d\n", DataArrays.nNode);

% Element order
DataArrays.ele_order = get_ele_order_1(msh);
fprintf("  Element order: %d\n", DataArrays.ele_order);

if DataArrays.ele_order == 1
    DataArrays.nLocBas = 3;
    DataArrays.nElem = msh.nbTriangles;
elseif DataArrays.ele_order == 2
    DataArrays.nLocBas = 6;
    DataArrays.nElem = msh.nbTriangles6;
end
fprintf("  Total elements: %d\n", DataArrays.nElem);

% Nodes' coordinate
fprintf("  Reading nodes ...\n");
DataArrays.Nodes = msh.POS';

% IEN array
fprintf("  Reading IEN ...\n");
DataArrays.IEN_v = get_IEN_v(msh, DataArrays.ele_order);

% Find Dirichlet nodes
fprintf("  Reading Dirichlet nodes...\n");
BCs.Dirichlet.DN = get_Dirichlet_node(msh, "Diri", DataArrays.ele_order);
BCs.Dirichlet.DiriFace = get_face2elem_faceid(msh, "Diri", DataArrays.IEN_v, DataArrays.ele_order);

% Construct ID array
fprintf("  Create ID array...\n");
DataArrays.ID = create_ID(DataArrays.nNode, DataArrays.dof);
% DataArrays.ID = modify_ID(DataArrays.ID, BCs.Dirichlet.DN, 1); 

% Construct LM array
fprintf("  Create LM array...\n");
DataArrays.LM = create_LM(DataArrays.IEN_v, DataArrays.ID);

fprintf("Preprocess finished.\n");

end

