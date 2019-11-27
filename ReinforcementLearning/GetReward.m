function reward=GetReward(S,cur_row,cur_col,rot_idx,action)
    % Computed the reward resulting from a selected action
    %
    % Input:
    % - S: n x m matrix that stores a current state
    % - cur_row: the row position of a robot's centroid
    % - cur_col: the column position of a robot's centroid
    % - rot_idx: a scalar indicating which rotation 
    % - action: a scalar indicating a selected action
    %      (1-going down, 2-going right, 3-going up, 4-going left
    %        5-rotating to rotation position 1, 6-rotating to rotation position 2
    %        7-rotating to rotation position 3, 8-rotating to rotation position 4).
    % Output:
    % - reward: a scalar indicating the received reward
   
   MIN_VAL=-99;
   MAX_VAL=99;

   % Getting next state
   [next_S,~,~,~]=MakeNextState(S,cur_row,cur_col,rot_idx,action);
   
   % Getting a reward
   if max(next_S(:))<=0 % CASE 1: the robot either hit the wall or stepped outside of the grid
       reward = MIN_VAL;
   else %the robot performed a valid action
      if isGoal(next_S) % CASE 2: case where the robot reached the goal
          reward = MAX_VAL;
      elseif  next_S == S% CASE 3: case where the next state is the same as the current state
          reward = MIN_VAL;
      else % CASE 4: every other case
          reward = 0;
      end
   end
   
end