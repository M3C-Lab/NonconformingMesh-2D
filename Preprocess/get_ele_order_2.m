function ele_order = get_ele_order_2(msh1, msh2)
% Get element order

if (msh1.nbTriangles > 0 && msh2.nbTriangles > 0)
    ele_order = 1;
elseif (msh1.nbTriangles6 > 0 && msh2.nbTriangles6 > 0)
    ele_order = 2;
else
    disp("Error, different element order of 2 domains");
    return;
end

end

