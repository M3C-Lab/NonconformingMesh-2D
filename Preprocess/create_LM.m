function LM = create_LM(IEN, ID)
% Create LM array

dof = size(ID, 1);
nLocBas = size(IEN, 1);
nElem = size(IEN, 2);

LM = zeros(dof * nLocBas, nElem);

for ee = 1 : nElem
    for aa = 1 : nLocBas
        node = IEN(aa, ee);
        for dd = 1 : dof
            LM((aa-1) * dof + dd, ee) = ID(dd, node);
        end
    end
end

end

