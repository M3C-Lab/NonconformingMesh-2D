function IEN_v = get_IEN_v(msh, ele_order)
% Get IEN array of volume element

if ele_order == 1
    IEN_v = msh.TRIANGLES(:, 1:3);
elseif ele_order == 2
    IEN_v = msh.TRIANGLES6(:, 1:6);
end
IEN_v = IEN_v';

end
