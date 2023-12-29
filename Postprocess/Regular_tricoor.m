function points = Regular_tricoor(N)
% To get some sampling points in a triangle with area coordinates .
% Input:
%   N: The parameter which controls the number of points in a triangle.
%       number of points = 0.5 * (N + 3) * (N + 2)
% Output:
%   points: The coordinates of the points.
%       points(:, ii): The area coordinates of the 'ii'th point.

points = [0.9999999999998, 0.0000000000001, 0.0000000000001;
          0.0000000000001, 0.9999999999998, 0.0000000000001;
          0.0000000000001, 0.0000000000001, 0.9999999999998];

if N > 0
    insert_points = zeros(3, 0.5 * (N + 3) * (N + 2) - 3);
    points = [points, insert_points];
    
    temp = 1;
    for ii = 1 : N
        p_left = points(:, 1) * (1 - ii / (N + 1)) + points(:, 2) * ii / (N + 1);
        p_right = points(:, 1) * (1 - ii / (N + 1)) + points(:, 3) * ii / (N + 1);
        
        for jj = 0 : ii
            insert_p = p_left * (1 - jj / ii) + p_right * jj / ii;
            points(:, temp + 3) = insert_p;
            temp = temp + 1;
        end
    end
    
    for jj = 1 : N
        insert_p = points(:, 2) * (1 - jj / (N + 1)) + points(:, 3) * (jj / (N + 1));
        points(:, temp + 3) = insert_p;
        temp = temp + 1;
    end
end

end

% EOF
