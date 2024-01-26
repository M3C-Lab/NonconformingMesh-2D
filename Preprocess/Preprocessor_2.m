function [DataArrays, BCs] = Preprocessor_2(msh1, msh2)
% Preprocessor for 2-domains
fprintf("\nPreprocess routine:\n")

% Degree of freedom at node
DataArrays.dof = 1;

% Number of nodes
DataArrays.num_node_1 = msh1.nbNod;
DataArrays.num_node_2 = msh2.nbNod;
DataArrays.nNode = msh1.nbNod + msh2.nbNod;
fprintf("  Nodes in mesh1: %d\n", msh1.nbNod);
fprintf("  Nodes in mesh2: %d\n", msh2.nbNod);
fprintf("  Total nodes: %d\n", DataArrays.nNode);

% Element order
DataArrays.ele_order = get_ele_order_2(msh1, msh2);
fprintf("  Element order: %d\n", DataArrays.ele_order);

if DataArrays.ele_order == 1
    DataArrays.nLocBas = 3;
    DataArrays.num_elem_1 = msh1.nbTriangles;
    DataArrays.num_elem_2 = msh2.nbTriangles;
elseif DataArrays.ele_order == 2
    DataArrays.nLocBas = 6;
    DataArrays.num_elem_1 = msh1.nbTriangles6;
    DataArrays.num_elem_2 = msh2.nbTriangles6;
end

fprintf("  Elements in mesh1: %d\n", DataArrays.num_elem_1);
fprintf("  Elements in mesh2: %d\n", DataArrays.num_elem_2);

% Number of elements
DataArrays.nElem = DataArrays.num_elem_1 + DataArrays.num_elem_2;
fprintf("  Total elements: %d\n", DataArrays.nElem);

% Nodes in domain_1
fprintf("  Reading nodes' coordinate of mesh1...\n");
nodes_1 = msh1.POS';

% IEN of volume element in domain_1
fprintf("  Reading IEN array of mesh1...\n");
DataArrays.IEN_v1 = get_IEN_v(msh1, DataArrays.ele_order);

% Nodes in domain_2
fprintf("  Reading nodes' coordinate of mesh2...\n");
nodes_2 = msh2.POS';

% IEN of volume element in domain_2
fprintf("  Reading IEN array of mesh2...\n");
DataArrays.IEN_v2 = get_IEN_v(msh2, DataArrays.ele_order);

% Find face2elem and face id
fprintf("  Reading face2elem and face id of interface1...\n");
BCs.traction.interface1 = get_face2elem_faceid(msh1, "Inter", DataArrays.IEN_v1, DataArrays.ele_order);
fprintf("    %d line elements on interface1.\n", size(BCs.traction.interface1, 2));

fprintf("  Reading face2elem and face id of interface2...\n");
BCs.traction.interface2 = get_face2elem_faceid(msh2, "Inter", DataArrays.IEN_v2, DataArrays.ele_order);
for ee = 1 : size(BCs.traction.interface2, 2)
    BCs.traction.interface2(1, ee) = BCs.traction.interface2(1, ee) + DataArrays.num_elem_1;
end
fprintf("    %d line elements on interface2.\n", size(BCs.traction.interface2, 2));

% Find Dirichlet nodes
fprintf("  Reading Dirichlet nodes...\n");
DiriNode1 = get_Dirichlet_node(msh1, "Diri", DataArrays.ele_order);
DiriFace1 = get_face2elem_faceid(msh1, "Diri", DataArrays.IEN_v1, DataArrays.ele_order);
fprintf("    %d nodes on Diri1.\n", length(DiriNode1));

DiriNode2 = get_Dirichlet_node(msh2, "Diri", DataArrays.ele_order);
for nn = 1 : length(DiriNode2)
    DiriNode2(nn) = DiriNode2(nn) + DataArrays.num_node_1;
end
DiriFace2 = get_face2elem_faceid(msh2, "Diri", DataArrays.IEN_v2, DataArrays.ele_order);
for ee = 1 : size(DiriFace2, 2)
    DiriFace2(1, ee) = DiriFace2(1, ee) + DataArrays.num_elem_1;
end
fprintf("    %d nodes on Diri2.\n", length(DiriNode2));
fprintf("  Combining Dirichlet node list ...\n");
BCs.Dirichlet.DN = [DiriNode1, DiriNode2];
BCs.Dirichlet.DiriFace = [DiriFace1, DiriFace2];

% All nodes
fprintf("  Combining nodes' coordinate...\n");
DataArrays.Nodes = [nodes_1, nodes_2];

% All elements
fprintf("  Combining IEN arrays...\n");
for ee = 1 : DataArrays.num_elem_2
    for aa = 1 : DataArrays.nLocBas
        DataArrays.IEN_v2(aa, ee) = DataArrays.IEN_v2(aa, ee) + DataArrays.num_node_1;
    end
end
DataArrays.IEN_v = [DataArrays.IEN_v1, DataArrays.IEN_v2];

% Construct ID array
fprintf("  Create ID array...\n");
DataArrays.ID = create_ID(DataArrays.nNode, DataArrays.dof);
% DataArrays.ID = modify_ID(DataArrays.ID, BCs.Dirichlet.DN, 1);

% Construct LM array
fprintf("  Create LM array...\n");
DataArrays.LM = create_LM(DataArrays.IEN_v, DataArrays.ID);

fprintf("Preprocess finished.\n");

end

