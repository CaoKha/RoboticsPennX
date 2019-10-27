function flag = triangle_intersection(P1, P2)
% triangle_test : returns true if the triangles overlap and false otherwise

% P1, P2: a 3 by 2 array (each), describing the vertices of a triangle,
% the first column corresponds to the x coordinates while the second column
% corresponds to the y coordinates

    function flag = sameside(p1,p2,a,b)
        % sameside : returns true if the p1,p1 lie on same sides of the
        % edge ab and false otherwise

        p1(3) = 0; p2(3) = 0; a(3) = 0; b(3) = 0;
        cp1 = cross(b-a, p1-a);
        cp2 = cross(b-a, p2-a);
        if(dot(cp1, cp2) >= 0)
            flag = true;
        else
            flag = false;
        end
    end

% Repeat the vertices for the loop
P1(4:5,:) = P1(1:2,:);
P2(4:5,:) = P2(1:2,:);

flag = true;

% Testing all the edges of P1
for i=1:3
    if(~sameside(P1(i,:), P2(1,:), P1(i+1,:), P1(i+2,:)) ...
            && sameside(P2(1,:), P2(2,:), P1(i+1,:), P1(i+2,:)) ...
            && sameside(P2(2,:), P2(3,:), P1(i+1,:), P1(i+2,:)))
        flag = false; return;
    end
end

% Testing all the edges of P2
for i=1:3
    if(~sameside(P2(i,:), P1(1,:), P2(i+1,:), P2(i+2,:)) ...
            && sameside(P1(1,:), P1(2,:), P2(i+1,:), P2(i+2,:)) ...
            && sameside(P1(2,:), P1(3,:), P2(i+1,:), P2(i+2,:)))
        flag = false; return;
    end
end

end