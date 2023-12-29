function sp_result = Sampler(fork, DataArrays, uh, n_sp, exact_u)
% To sample some points in each element, and get the result data of them.
fprintf("  Sampling:\n");
fprintf("    %d sampling points in each element ...\n", n_sp)
sampling_points = Regular_tricoor(n_sp);

if fork == 0
    fprintf("    Sampling in the whole domain ...\n");
    nElem = DataArrays.nElem;
    IEN_v = DataArrays.IEN_v;
elseif fork == 1
    fprintf("    Sampling in the left domain ...\n");
    nElem = DataArrays.num_elem_1;
    IEN_v = DataArrays.IEN_v1;
elseif fork == 2
    fprintf("    Sampling in the right domain ...\n");
    nElem = DataArrays.num_elem_2;
    IEN_v = DataArrays.IEN_v2;
else
    disp("Error, wrong fork.");
    return;
end

nb_sp = size(sampling_points, 2);
total_sp = nElem * nb_sp;

sp_result.Location = zeros(total_sp, 3);
sp_result.u_interpolated = zeros(total_sp, 1);

% Integral to calculate relative error
sp_result.Int_sq_u = 0.0;
sp_result.Int_sq_u_minus_uh = 0.0;

sp_result.h_meshsize = 0;
% The max mesh size.
node_ele = zeros(3, DataArrays.nLocBas);
uh_ele = zeros(DataArrays.dof * DataArrays.nLocBas, 1);
for ee = 1 : nElem
    for aa = 1 : DataArrays.nLocBas
        node_id = IEN_v(aa, ee);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
        uh_ele(aa) = uh(node_id);
    end
    
    h_ele = 2 * circumcircle([node_ele(1:2, 1), node_ele(1:2, 2), node_ele(1:2, 3)], 0);
    if h_ele > sp_result.h_meshsize
        sp_result.h_meshsize = h_ele;
    end
    
    for ss = 1 : nb_sp
        sp = sampling_points(:, ss); % The parametric coordinate of a sampling point.
        
        triquaddata = TriBasis(DataArrays.ele_order, node_ele, sp);

        sp_xyz = node_ele * triquaddata.Basis;
        % The coordinates of sampling points
        sp_result.Location(nb_sp * (ee - 1) + ss, 1) = sp_xyz(1);
        sp_result.Location(nb_sp * (ee - 1) + ss, 2) = sp_xyz(2);
        sp_result.Location(nb_sp * (ee - 1) + ss, 3) = ee;

        sp_result.u_interpolated(nb_sp * (ee - 1) + ss) = dot(uh_ele, triquaddata.Basis);
    end

    for qua = 1 : DataArrays.postquad.nqp
        triquaddata = TriBasis(DataArrays.ele_order, node_ele, DataArrays.postquad.qp(:, qua));

        qua_xyz = node_ele * triquaddata.Basis;

        u_qua = exact_u(qua_xyz(1), qua_xyz(2));
        uh_qua = dot(uh_ele, triquaddata.Basis);

        sp_result.Int_sq_u = sp_result.Int_sq_u + DataArrays.postquad.wq(qua) * norm(triquaddata.n_out) * u_qua^2;
        sp_result.Int_sq_u_minus_uh = sp_result.Int_sq_u_minus_uh + DataArrays.postquad.wq(qua) * norm(triquaddata.n_out)* (u_qua - uh_qua)^2;
    end
end

fprintf("  Sampling finished.\n");
return;
end

% EOF

