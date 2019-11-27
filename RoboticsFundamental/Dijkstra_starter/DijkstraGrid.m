function [route,numExpanded] = DijkstraGrid(input_map, start_coords, dest_coords, drawMap)
% Run Dijkstra's algorithm on a grid.
% Inputs : 
%   input_map : a logical array where the freespace cells are false or 0 and
%   the obstacles are true or 1
%   start_coords and dest_coords : Coordinates of the start and end cell
%   respectively, the first entry is the row and the second the column.
% Output :
%    route : An array containing the linear indices of the cells along the
%    shortest route from start to dest or an empty array if there is no
%    route. This is a single dimensional vector
%    numExpanded: Remember to also return the total number of nodes
%    expanded during your search. Do not count the goal node as an expanded node.


% set up color map for display
% 1 - white - clear cell
% 2 - black - obstacle
% 3 - red = visited
% 4 - blue  - on list
% 5 - green - start
% 6 - yellow - destination

cmap = [1 1 1; ...
        0 0 0; ...
        1 0 0; ...
        0 0 1; ...
        0 1 0; ...
        1 1 0; ...
	0.5 0.5 0.5];

colormap(cmap);

% variable to control if the map is being visualized on every
% iteration
drawMapEveryTime = drawMap;

[nrows, ncols] = size(input_map);

% map - a table that keeps track of the state of each grid cell
map = zeros(nrows, ncols);

map(~input_map) = 1;   % Mark free cells
map(input_map)  = 2;   % Mark obstacle cells

% Generate linear indices of start and dest nodes
start_node = sub2ind(size(map), start_coords(1), start_coords(2));
dest_node  = sub2ind(size(map), dest_coords(1),  dest_coords(2));

map(start_node) = 5;
map(dest_node)  = 6;

% Initialize distance array
distanceFromStart = Inf(nrows, ncols);

% For each grid cell this array holds the index of its parent
parent = zeros(nrows, ncols);

distanceFromStart(start_node) = 0;

% keep track of number of nodes expanded 
numExpanded = 0;

% Main Loop
while true
    
    % Draw current map
    map(start_node) = 5;
    map(dest_node) = 6;
    
    % make drawMapEveryTime = true if you want to see how the 
    % nodes are expanded on the grid. 
    if (drawMapEveryTime)
        image(1.5, 1.5, map);
        grid on;
        axis image;
        drawnow;
    end
    
    % Find the node with the minimum distance
    [min_dist, current] = min(distanceFromStart(:));
    if ((current == dest_node) || isinf(min_dist))
        break
    end       
 % Update map
    map(current) = 3;         % mark current node as visited  
    
          
    % Compute row, column coordinates of current node
    [i, j] = ind2sub(size(distanceFromStart), current);

    
    
   % ********************************************************************* 
    % YOUR CODE BETWEEN THESE LINES OF STARS
    % find neighbors
    neigh(1:4,1:2) = [i+[-1;0;0;+1] j+[0;-1;+1;0] ];
    % remove all neighbors that have index 0
    sz = size(distanceFromStart);
    neigh = neigh(all(neigh,2) & neigh(:,1) <= sz(1) & neigh(:,2) <= sz(2),:);
    % check if neighbor nodes were visited and check their distance from start,
    % if it = Inf, update new distance = distance of current node + 1.
    % after that, update the node as visited
    nei_r = size(neigh,1);
    for di = 1:nei_r
        if map(neigh(di,1),neigh(di,2)) ~= 2 && map(neigh(di,1),neigh(di,2)) ~= 5 && map(neigh(di,1),neigh(di,2)) ~= 3
            distance_d = distanceFromStart(i,j) +  sqrt((neigh(di,1) - i)^2 + (neigh(di,2)-j)^2);
            if distance_d < distanceFromStart(neigh(di,1),neigh(di,2))
                distanceFromStart(neigh(di,1),neigh(di,2)) = distance_d;
                parent(neigh(di,1),neigh(di,2)) = current;  
            end
        end
    end
    distanceFromStart(current) = Inf; % remove this node from further consideration 
 
           
            
            
            
    % Visit each neighbor of the current node and update the map, distances
    % and parent tables appropriately.
    
    

    %*********************************************************************
end
for row = 1:10
    for col = 1:10
        if ((start_coords(1) - dest_coords(1))^2 + (start_coords(2) - dest_coords(2))^2) ~= 0
            if map(row,col) == 3 || map(row,col) == 5
                numExpanded = numExpanded + 1;
            end
        else 
            numExpanded = 1;
        end
    end
end
%% Construct route from start to dest by following the parent links
if (isinf(distanceFromStart(dest_node)))
    route = [];
else
    route = dest_node;
    
    while (parent(route(1)) ~= 0)
        route = [parent(route(1)), route];
    end
    
    % Snippet of code used to visualize the map and the path
    for k = 2:length(route) - 1        
        map(route(k)) = 7;
        pause(0.1);
        image(1.5, 1.5, map);
        grid on;
        axis image;
    end
end

end
