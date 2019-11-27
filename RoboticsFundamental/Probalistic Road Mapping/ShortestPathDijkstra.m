function route = ShortestPathDijkstra (edges, edge_lengths, start, dest)
% This function finds the shortest path through a graph using Dijkstra's
% algorithm. The first input, map is a sparse matrix whose entries
% indicate the distance values all of which must be greater than zero
% source and dest indicate the indices of the start and destination nodes.


% Construct a sparse matrix using the information in edges and
% edge_distances

if (any(edge_lengths <= 0))
    error ('Edge Distances cannot be less than or equal to zero');
end

i = edges(:,1);
j = edges(:,2);

map = sparse([i,j], [j,i], [edge_lengths, edge_lengths]);

% Figure out the number of nodes
n = size(map,1);

if ((start > n) || (dest > n))
    % This could happen if the start or end node has no edges
    route = [];
    return;
end

% Initialization
distances = Inf(n,1);
parent = NaN(n,1);

distances(start) = 0;
parent(start) = 0;

while true
    % Find the node with the smallest distance
    [min_dist,idx] = min(distances);
    
    if ((idx == dest) || isinf(min_dist))
       break;
    end
    
    i = find(map(idx,:));
    
    % ********************************************************************* 
    % YOUR CODE BETWEEN THESE LINES OF STARS
    sizeofneighbor = size(i,2);
    for i_index = 1:sizeofneighbor
        if isnan(distances(i(i_index))) == false 
            alt = distances(idx) + map(idx,i(i_index));
            if alt < distances(i(i_index))
                distances(i(i_index)) = alt;
                parent(i(i_index)) = idx;
            end
        end
    end
    distances(idx) = NaN; % remove this node from further consideration 
    %*********************************************************************

end

% Construct the route
if (isinf(distances(dest)))
    route = [];
else
    route = [dest];
    
    while (parent(route(1)) ~= 0)
        route = [parent(route(1)), route];
    end
end

