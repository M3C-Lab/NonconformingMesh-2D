function [K_loc_loc, K_oppo_oppo, K_loc_oppo, K_oppo_loc] = LocAssem_NitscheInterface(local_basis, n1, oppo_basis, DataArrays, Paras)

local_dof = DataArrays.nLocBas * DataArrays.dof;
K_loc_loc = zeros(local_dof, local_dof);
K_oppo_oppo = zeros(local_dof, local_dof);
K_loc_oppo = zeros(local_dof, local_dof);
K_oppo_loc = zeros(local_dof, local_dof);

for A = 1 : DataArrays.nLocBas
    NA1 = local_basis.Basis(A);
    NA1_x = local_basis.dN_dx(A);
    NA1_y = local_basis.dN_dy(A);

    NA2 = oppo_basis.Basis(A);
    NA2_x = oppo_basis.dN_dx(A);
    NA2_y = oppo_basis.dN_dy(A);

    for B = 1 : DataArrays.nLocBas
        NB1 = local_basis.Basis(B);
        NB1_x = local_basis.dN_dx(B);
        NB1_y = local_basis.dN_dy(B);
    
        NB2 = oppo_basis.Basis(B);
        NB2_x = oppo_basis.dN_dx(B);
        NB2_y = oppo_basis.dN_dy(B);

        K_loc_loc(A, B) = K_loc_loc(A, B) ...
            - Paras.adjoint_coef * (NA1 * (NB1_x * n1(1) + NB1_y * n1(2)) + (NA1_x * n1(1) + NA1_y * n1(2)) * NB1) ...
            + Paras.penalty_coef * NA1 * NB1;

        K_oppo_oppo(A, B) = K_oppo_oppo(A, B) ...
            - Paras.adjoint_coef * (NA2 * (NB2_x * n1(1) + NB2_y * n1(2)) + (NA2_x * n1(1) + NA2_y * n1(2)) * NB2) ...
            + Paras.penalty_coef * NA2 * NB2;

        K_loc_oppo(A, B) = K_loc_oppo(A, B) ...
            + Paras.adjoint_coef * (NA1 * (NB2_x * n1(1) + NB2_y * n1(2)) + (NA1_x * n1(1) + NA1_y * n1(2)) * NB2) ...
            - Paras.penalty_coef * NA1 * NB2;

        K_oppo_loc(A, B) = K_oppo_loc(A, B) ...
            + Paras.adjoint_coef * (NA2 * (NB1_x * n1(1) + NB1_y * n1(2)) + (NA2_x * n1(1) + NA2_y * n1(2)) * NB1) ...
            - Paras.penalty_coef * NA2 * NB1;

    end
end

end

