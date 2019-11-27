function is_goal=isGoal(S)
    % Check if a current state is a goal state
    %
    % Input:
    % - S: n x m matrix that stores a current state
    % Output:
    % - is_goal: a binary variable with a value of 1 if S is the goal state, and 0 otherwise
    
    if S == MakeState(7,7,[2 4; 3 4; 4 4; 5 4],6,7,3)
        is_goal=1;
    else
        is_goal = 0;
    end
end