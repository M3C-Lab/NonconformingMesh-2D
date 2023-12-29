function Initial_ID = create_ID(num_node, dof)
% Create an initial ID array
temp = 1;
Initial_ID = zeros(dof, num_node);
for nn = 1 : num_node
    for ii = 1 : dof
        Initial_ID(ii, nn) = temp;
        temp = temp + 1;
    end
end

end

