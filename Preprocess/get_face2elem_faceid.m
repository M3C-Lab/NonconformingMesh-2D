function volid_faceid = get_face2elem_faceid(msh, surface_name, IEN_v, ele_order)
% Get the vol_ele id and face_id of each surface element

% Number of physical groups
num_PG = size(msh.PhyGrp, 1);

% Find surface id
for ii = 1 : num_PG
    if strcmp(surface_name, msh.PhyGrp{ii, 3})
        PG_id = msh.PhyGrp{ii, 2};
        break;
    end
    if ii == num_PG
        disp("Error, wrong surface name, cannot find physical group.");
        return;
    end
end

if ele_order == 1
    num_total_sur_ele = msh.nbLines;
    sur_nLocBas = 2;
    sur_ele = msh.LINES;
elseif ele_order == 2
    num_total_sur_ele = msh.nbLines3;
    sur_nLocBas = 3;
    sur_ele = msh.LINES3;
end

% Find surface IEN
IEN_s = [];
for ee = 1 : num_total_sur_ele
    if sur_ele(ee, sur_nLocBas + 1) == PG_id    % Catch surface element
        localien = sur_ele(ee, 1 : sur_nLocBas);
        localien = localien';
        IEN_s = [IEN_s, localien];  % Insert local IEN
    end
end

num_sur_ele = size(IEN_s, 2);
num_vol_ele = size(IEN_v, 2);

volid_faceid = zeros(2, num_sur_ele);
for ee = 1 : num_sur_ele
    sur_vertices = IEN_s(1:2, ee);
    for EEE = 1 : num_vol_ele
        vol_vertices = IEN_v(1:3, EEE);
        for mm = 1 : 3
            if vol_vertices(mm) == sur_vertices(1)
                for nn = 1 : 3
                    if vol_vertices(nn) == sur_vertices(2)
                        volid_faceid(1, ee) = EEE;
                        volid_faceid(2, ee) = 6 - mm - nn;
                    end
                    if volid_faceid(1, ee)~= 0
                        break;
                    end
                end
            end
            if volid_faceid(1, ee)~= 0
                break;
            end
        end
    end
end

end

