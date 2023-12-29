function ele_order = get_ele_order_1(msh)
% Get element order for 1-domain

if (msh.nbTriangles > 0)
    ele_order = 1;
elseif (msh.nbTriangles6 > 0)
    ele_order = 2;
else
    disp("Error, unknown element.");
    return;
end

end

