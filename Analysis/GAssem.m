 function [K, F] = GAssem(K, F, DataArrays, f)
% Main global assembly

fprintf("  Global assembling in the whole domain ...\n");
for ee = 1 : DataArrays.nElem
    [K_ele, F_ele] = LocAssem(ee, DataArrays, f);

    local_dof = size(DataArrays.LM, 1);
    for aa = 1 : local_dof
        LM_a = DataArrays.LM(aa, ee);
        if LM_a > 0
            F(LM_a) = F(LM_a) + F_ele(aa);
            for bb = 1 : local_dof
                LM_b = DataArrays.LM(bb, ee);
                if LM_b > 0
                    K(LM_a, LM_b) = K(LM_a, LM_b) + K_ele(aa, bb);
                    % else do nothing since Dirichlet BC g=0
                end
            end
        end
    end
end

end

