function action=pickAction(S,cur_row,cur_col,rot_idx,nn,epsilon)
    % Pick an action for the robot to carry out
    %
    % Input:
    % - S: n x m matrix that stores a current state
    % - cur_row: the row position of a robot's centroid
    % - cur_col: the column position of a robot's centroid
    % - rot_idx: a scalar indicating which rotation 
    % - nn: a variable storing a neural network
    % - epsilon: a parameter that controls how often we select action randomly
    % Output:
    % - action: a scalar value in the interval [1,8] (inclusive) that depicts which action we selected
    
    rand('state',0);
    
    %% Picking Action
    rand_val=rand();
    
    % Feed State S to neural network to receive probability of 8 actions at
    % the last layer
    nn = ForwardPass(nn,S(:)');
    
    if rand_val>epsilon
        % Use Neural Network to select an action
        [~,action] = max(nn.a{end});     
    else
        % Select an action randomly
        action= randi([1 8]);
    end
   
end
