function nn=DeepQLearning(action,reward,S,new_S,nn)
    % Update the weights of the network
    %
    % Input:
    % - action: a scalar indicating which action was selected in state S
    % - reward: a scalar indicating the reward value that was obtained by transition to state new_S
    % - S: n x m matrix that stores the previous state
    % - new_S: n x m matrix that stores the current state
    % - nn: a variable storing a neural network
    % Output:
    % - nn: a variable storing a neural network with the updated weights after a single Q-learning step.
    
    
    lr=0.01; %learning rate
    gamma=0.8;
    n_layers=numel(nn.W)+1;
    
    %% Initializing target variable
    nn = ForwardPass(nn, S(:)');
    target=nn.a{n_layers}';

    %% Computing the maximum Q value of next action by 
    %% Feeding the new state new_S through a network
    nn = ForwardPass(nn,new_S(:)');
    
    %% Setting a target(action) value
    target(action)= reward + gamma*max(nn.a{n_layers});   
    
    %% Performing the Backward pass
    nn = ForwardPass(nn, S(:)');
    nn = BackwardPass(nn,target);
    
    %% Updating the weights
    for l = 1 : (n_layers - 1)
       nn.W{l} = nn.W{l} - lr*nn.gradW{l};   
    end

end