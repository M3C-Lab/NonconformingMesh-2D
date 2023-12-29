addpath("Preprocess/", "Analysis/", "Element_and_Quadrature/", "Postprocess/");

clc; clear all; close all;
% exact solution: u = x * (1 - x) * y * (1 - y)
u = @(x, y)  x * (1 - x) * y * (1 - y);

% exact "f = -Laplacian u"
f = @(x, y)  2 * (y * (1 - y) + x * (1 - x));

%% Preprocess
msh = load_gmsh2("Mesh/1domain.msh");
[DataArrays, BCs] = Preprocessor_1(msh);

% Set quadrature rules
if DataArrays.ele_order == 1
    [DataArrays.triquad.qp, DataArrays.triquad.wq, DataArrays.triquad.nqp] = TriQuad(3);
elseif DataArrays.ele_order == 2
    [DataArrays.triquad.qp, DataArrays.triquad.wq, DataArrays.triquad.nqp] = TriQuad(6);
end
[DataArrays.postquad.qp, DataArrays.postquad.wq, DataArrays.postquad.nqp] = TriQuad(6);

%% Analysis
% Total degree of freedom
total_dof = DataArrays.dof * DataArrays.nNode;

% Initilize K and F
K = sparse(total_dof, total_dof);
F = zeros(total_dof, 1);

fprintf("\nGlobal assembly routine:\n");

% Global assembly
% Main global assembly
[K, F] = GAssem(K, F, DataArrays, f);

fprintf("Global assembly finifshed.\n");

fprintf("\nApply essential BC ...\n");
[K, F] = EssBC(K, F, total_dof, BCs);

fprintf("Solving ...\n");
uh = K \ F;

%% Postprocess
fprintf("\nPostprocess routine:\n");
sp_result = Sampler(0, DataArrays, uh, 5, u);

fprintf("  Plotting ...\n")
PlotSampling(sp_result, "Numerical solution", 1);
PlotMesh(DataArrays, sp_result.h_meshsize, 2, 0);

hh = sp_result.h_meshsize;
rel_error = sp_result.Int_sq_u_minus_uh / sp_result.Int_sq_u;
rel_error = sqrt(rel_error);

fprintf("Done!\n");
