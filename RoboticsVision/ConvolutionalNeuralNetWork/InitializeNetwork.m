% architecture=[336 100 20];
% nn = InitializeNetwork(architecture);

function nn = InitializeNetwork(architecture)
    % initialize the network
    %
    % Input:
    % - architecture: an 1 x l dimensional vector indicating a number of nodes at each layer of a network. For instance, architecture(1) is the feature dimensionality, architecture(2) is the number of nodes in the first hidden layer, and so on.
    % Output:
    % - nn: a neural network structure storing randomly initialized parameters of a network in the variables nn.W{1},nn.W{2},etc.
    
    rand('state',0)  
    for i = 2 : size(architecture,2)
        rows = architecture(i); %specificy the number of rows for the parameter matrix;
        cols = architecture(i-1); %specify the number of columns in the parameter matrix;
        nn.W{i-1}= (rand(rows, cols) - 0.5) * 2 * 4 * sqrt(6 / (rows + cols));
    end
end