function nn = BackwardPass(nn,y)
    % perform the forward pass inside the network
    %
    % Input:
    % - nn: a structure storing the parameters of the network
    % - y: a ground truth indicator matrix, where y(i,j)=1 indicates that a data point i belongs to an object class j.
    % Output:
    % - nn: a new neural network variable where the values nn.gradZ{l} and nn.gradW{l} are updated for every layer of the network.
    
   
    %number of layers in a network
    n_layers=numel(nn.W)+1;
    
    %computing the overall error of the network
    error= nn.a{n_layers} - y';
    
    %computing the error gradient in the last layer of the network
    nn.gradZ{n_layers} = error .* (nn.a{n_layers} .* (1 - nn.a{n_layers}));
    
    %looping through the layers backwards
    for i = (n_layers-1) : -1 : 2
        
        % 1) Derivative of the Sigmoid function                     gradSigm(l) = da(l)/dz(l)
        nn.gradSigm{i} = nn.a{i}.*(1-nn.a{i});
        
        % 2) Derivative with respect to activation units a{i}       gradA(l) = dL(l)/da(l+1)
        nn.gradA{i} = (nn.W{i}.')*nn.gradZ{i+1};
        
        % 3) Derivative with respect to activation units z{i}       gradZ(l) = dL(l)/dz(l)           
        nn.gradZ{i} = nn.gradA{i}.*nn.gradSigm{i};  

    end

    % Computing the gradients with respect to the fully connected layer parameters
    for i = 1 : (n_layers - 1)
        
        % 4) Derivative with respect to the fully connected layer parameters W{i}
        nn.gradW{i} = nn.gradZ{i+1}*(nn.a{i}.');
        
        % Averaging over all examples in the batch
        nn.gradW{i}=nn.gradW{i}/size(y,1);
        
    end
end