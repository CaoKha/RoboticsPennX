% %% Create a starting state
% rows=7; cols=7;
% walls=[2 4; 3 4; 4 4; 5 4];
% cur_row=2; cur_col=1; rot_idx=1;
% start_S=MakeState(rows,cols,walls,cur_row,cur_col,rot_idx);

% % Current state
% rows=7; cols=7;
% walls=[2 4; 3 4; 4 4; 5 4];
% cur_row=2; cur_col=1; rot_idx=1;
% S=MakeState(rows,cols,walls,cur_row,cur_col,rot_idx);
% 
% % Neural Network
% nn = InitializeNetwork([rows*cols 8]);
% 
% % Picking an action
% epsilon=0.5;
% action=pickAction(S,cur_row,cur_col,rot_idx,nn,epsilon);

% % Current state
% rows=7; cols=7;
% walls=[2 4; 3 4; 4 4; 5 4];
% cur_row=2; cur_col=1; rot_idx=1;
% S=MakeState(rows,cols,walls,cur_row,cur_col,rot_idx);
% 
% % Neural Network
% nn = InitializeNetwork([rows*cols 8]);
% 
% % Picking action deterministically
% epsilon=0;
% action=pickAction(S,cur_row,cur_col,rot_idx,nn,epsilon);
% 
% % Getting a reward
% reward=GetReward(S,cur_row,cur_col,rot_idx,action);
% 
% %Transitioning to a New State
% [new_S,~,~,~]=MakeNextState(S,cur_row,cur_col,rot_idx,action);
% 
% % Deep Q-Learning
% nn=DeepQLearning(action,reward,S,new_S,nn);

nn=ReinforcementLearning();
action_path=Prediction(nn);
number_of_actions_taken=size(action_path,1);