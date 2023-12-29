function BCs = get_seeking_priority(DataArrays, BCs)
% Get the seeking priority for interface element

% Number of elements on interface1
num_inter_elem1 = size(BCs.traction.interface1, 2);

centroid_interface1 = zeros(3, num_inter_elem1);

for ee = 1 : num_inter_elem1
    vol_ele_id = BCs.traction.interface1(1, ee);
    line_id = BCs.traction.interface1(2, ee);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end

    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, DataArrays.linequad.qp);

    centroid_interface1(:, ee) = get_line_centroid(node_line, DataArrays);
end

% Number of elements on interface2
num_inter_elem2 = size(BCs.traction.interface2, 2);

centroid_interface2 = zeros(3, num_inter_elem2);

for ee = 1 : num_inter_elem2
    vol_ele_id = BCs.traction.interface2(1, ee);
    line_id = BCs.traction.interface2(2, ee);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end

    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, DataArrays.linequad.qp);

    centroid_interface2(:, ee) = get_line_centroid(node_line, DataArrays);
end

BCs.seeking_priority1 = zeros(num_inter_elem2, num_inter_elem1);
for ee = 1 : num_inter_elem1
    centroid_here = centroid_interface1(:, ee);
    priority_list = [1];
    for EEE = 2 : num_inter_elem2
        dist = norm(centroid_interface2(:, EEE) - centroid_here);
        
        currLength = length(priority_list);

        location = currLength + 1;

        % Find a location to insert this opposite element
        for comp = 1 : currLength
            dist_comp = norm(centroid_interface2(:, priority_list(comp)) - centroid_here);
            if dist < dist_comp
                location = comp;
                break;
            end
        end

        if location == 1
            priority_list = [EEE; priority_list];
        elseif location == currLength + 1
            priority_list = [priority_list; EEE];
        else
            priority_list = [priority_list(1 : location - 1, 1); EEE; priority_list(location : currLength, 1)];
        end
    end

    BCs.seeking_priority1(:, ee) = priority_list;
end

BCs.seeking_priority2 = zeros(num_inter_elem1, num_inter_elem2);
for ee = 1 : num_inter_elem2
    centroid_here = centroid_interface2(:, ee);
    priority_list = [1];
    for EEE = 2 : num_inter_elem1
        dist = norm(centroid_interface1(:, EEE) - centroid_here);

        currLength = length(priority_list);

        location = currLength + 1;

        % Find a location to insert this opposite element
        for comp = 1 : currLength
            dist_comp = norm(centroid_interface1(:, priority_list(comp)) - centroid_here);
            if dist < dist_comp
                location = comp;
                break;
            end
        end

        if location == 1
            priority_list = [EEE; priority_list];
        elseif location == currLength + 1
            priority_list = [priority_list; EEE];
        else
            priority_list = [priority_list(1 : location - 1, 1); EEE; priority_list(location : currLength, 1)];
        end
    end

    BCs.seeking_priority2(:, ee) = priority_list;
end

end

