function centroid = get_line_centroid(node_line, DataArrays)
% Get the centroid of a line element

x_Length = 0.0; y_Length = 0.0; z_Length = 0.0; Length = 0.0;

for qua = 1 : DataArrays.linequad.nqp
    line = LineBasis(DataArrays.ele_order, node_line, DataArrays.linequad.qp(qua));

    qua_coor = node_line * line.Basis;

    J = norm(line.tangent_vector);

    x_Length = x_Length + DataArrays.linequad.wq(qua) * qua_coor(1) * J;
    y_Length = y_Length + DataArrays.linequad.wq(qua) * qua_coor(2) * J;
    z_Length = z_Length + DataArrays.linequad.wq(qua) * qua_coor(3) * J;

    Length = Length + DataArrays.linequad.wq(qua) * J;
end

centroid = [x_Length; y_Length; z_Length] / Length;

end

