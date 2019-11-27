function S=MakeState(n,m,walls,cur_row,cur_col,rot_idx)
    % Create a current state
    %
    % Input:
    % - n: number of rows in the environment grid
    % - m: number of columns in the environment grid
    % - walls: d x 2 matrix containing (row,column) positions of the walls in the grid
    % - cur_row: the row position of a robot's centroid
    % - cur_col: the column position of a robot's centroid
    % - rot_idx: a scalar indicating which rotation position the robot is currently in
    % Output:
    % - S: n x m matrix storing the state: wall locations should have the values of -1, every position occupied by a robot should have values
    % of 1, the rest of the gridd cells should have zero values. If the
    % provided state is not valid then all the values in S should have zero
    % values
    
    S=zeros(n,m);  
    
    % Fill in The Wall values
    for i = 1:size(walls,1)
        S(walls(i,1),walls(i,2)) = -1;
    end
    
    % Find the locations in the grid occupied by a robot
    robot_idx = zeros(5,2);
    if rot_idx==1
        robot_idx(1,1) = cur_row-1;
        robot_idx(1,2) = cur_col;
        robot_idx(2,1) = cur_row-1;
        robot_idx(2,2) = cur_col+1;
        robot_idx(3,1) = cur_row;
        robot_idx(3,2) = cur_col;
        robot_idx(4,1) = cur_row+1;
        robot_idx(4,2) = cur_col;
        robot_idx(5,1) = cur_row+1;
        robot_idx(5,2) = cur_col+1;
                            
    elseif rot_idx==2
        robot_idx(1,1) = cur_row;
        robot_idx(1,2) = cur_col-1;
        robot_idx(2,1) = cur_row;
        robot_idx(2,2) = cur_col;
        robot_idx(3,1) = cur_row;
        robot_idx(3,2) = cur_col+1;
        robot_idx(4,1) = cur_row+1;
        robot_idx(4,2) = cur_col-1;
        robot_idx(5,1) = cur_row+1;
        robot_idx(5,2) = cur_col+1;

    elseif rot_idx==3
        robot_idx(1,1) = cur_row-1;
        robot_idx(1,2) = cur_col-1;
        robot_idx(2,1) = cur_row-1;
        robot_idx(2,2) = cur_col;
        robot_idx(3,1) = cur_row;
        robot_idx(3,2) = cur_col;
        robot_idx(4,1) = cur_row+1;
        robot_idx(4,2) = cur_col-1;
        robot_idx(5,1) = cur_row+1;
        robot_idx(5,2) = cur_col;

    elseif rot_idx==4
        robot_idx(1,1) = cur_row-1;
        robot_idx(1,2) = cur_col-1;
        robot_idx(2,1) = cur_row-1;
        robot_idx(2,2) = cur_col+1;
        robot_idx(3,1) = cur_row;
        robot_idx(3,2) = cur_col-1;
        robot_idx(4,1) = cur_row;
        robot_idx(4,2) = cur_col;
        robot_idx(5,1) = cur_row;
        robot_idx(5,2) = cur_col+1;
    end
    
    % Define a state variable to check if S is valid
    state = 1;
    % Check if the robot is fully within the boundaries of a grid
    for i = 1:size(robot_idx,1)
        if (robot_idx(i,1) <= 0) || (robot_idx(i,1) > n) || (robot_idx(i,2) <= 0) || (robot_idx(i,2) > m)
            state = 0;
            break;
        end
    end
    
    % Check if the robot's location doesn't overlap with the locations of the walls
    
    % Check if the robot is fully within the boundaries of a grid
    for i = 1:size(robot_idx,1)
        for j = 1:size(walls,1)
            if (robot_idx(i,1) == walls(j,1)) && (robot_idx(i,2) == walls(j,2))
                state = 0;
                break;
            end
        end
    end
    
    % If the state is not valid return S with all zero values, otherwise return the state S filled in with the correct values
    if state == 0
        S = zeros(n,m);
    else
        for i = 1:size(robot_idx,1)
            S(robot_idx(i,1),robot_idx(i,2)) = 1;
        end
    end
end