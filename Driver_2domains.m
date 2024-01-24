addpath("Preprocess/", "Analysis/", "Element_and_Quadrature/", "Postprocess/");

clc; clear all; close all;
% exact solution: u = x * (1 - x) * y * (1 - y)
u = @(x, y)  x * (1 - x) * y * (1 - y);

% exact "f = -Laplacian u"
f = @(x, y)  2 * (y * (1 - y) + x * (1 - x));

%% Preprocess
msh1 = load_gmsh2("Mesh/2domains_left.msh");
msh2 = load_gmsh2("Mesh/2domains_right.msh");

[DataArrays, BCs] = Preprocessor_2(msh1, msh2);

% Set quadrature rules
if DataArrays.ele_order == 1
    [DataArrays.triquad.qp, DataArrays.triquad.wq, DataArrays.triquad.nqp] = TriQuad(3);
    [DataArrays.linequad.qp, DataArrays.linequad.wq] = Gauss(3, 0.0, 1.0);
    DataArrays.linequad.nqp = 3;
elseif DataArrays.ele_order == 2
    [DataArrays.triquad.qp, DataArrays.triquad.wq, DataArrays.triquad.nqp] = TriQuad(6);
    [DataArrays.linequad.qp, DataArrays.linequad.wq] = Gauss(6, 0.0, 1.0);
    DataArrays.linequad.nqp = 6;
end
[DataArrays.postquad.qp, DataArrays.postquad.wq, DataArrays.postquad.nqp] = TriQuad(6);

% Set parameters for Nitsche's method
Paras.penalty_coef = 100;
Paras.adjoint_coef = 1.0;

% Set opposite seeking priority for interface elements
BCs = get_seeking_priority(DataArrays, BCs);

%% Analysis
% Total degree of freedom
total_dof = DataArrays.dof * DataArrays.nNode;

% Initilize K and F
K = sparse(total_dof, total_dof);
F = zeros(total_dof, 1);

fprintf("\nGlobal assembly routine:\n");

% Global assembly
% Main global assembly
[K, F] = GAssem(K, F, DataArrays, u, f);

% Global assembly for interface
[K, F] = GAssem_TryNitsche(K, F, DataArrays, BCs, Paras, u);

fprintf("Global assembly finifshed.\n");

fprintf("\nApply essential BC ...\n");
[K, F] = EssBC(K, F, total_dof, BCs);

fprintf("Solving ...\n");
uh = K \ F;

%% Postprocess
fprintf("\nPostprocess routine:\n");
sp_result1 = Sampler(1, DataArrays, uh, 5, u);
sp_result2 = Sampler(2, DataArrays, uh, 5, u);

fprintf("  Plotting ...\n")
PlotSampling(sp_result1, "Numerical solution in left domain", 1);
PlotMesh(DataArrays, sp_result1.h_meshsize, 2, 1);

PlotSampling(sp_result2, "Numerical solution in right domain", 3);
PlotMesh(DataArrays, sp_result2.h_meshsize, 4, 2);

hh = min(sp_result1.h_meshsize, sp_result2.h_meshsize);
rel_error = (sp_result1.Int_sq_u_minus_uh + sp_result2.Int_sq_u_minus_uh) / (sp_result1.Int_sq_u + sp_result2.Int_sq_u);
rel_error = sqrt(rel_error);

fprintf("Done!\n");