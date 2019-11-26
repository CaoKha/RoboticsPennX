function preds=PredictNetwork(nn,test_x)
    % produce the predictions of the trained network
    %
    % Input:
    % - nn: a structure storing the parameters of the network
    % - test_x: n x d dimensional feature matrix storing d dimensional feature for n data points
    % Output:
    % - preds: n x 1 matrix that stores the predicted object class indices ranging from values [1...20].

    %1) Perform the forward pass
    nn = ForwardPass(nn,test_x);
    %2) Select object classes cororesponding to maximum probabilities and store them into 'preds' variable
    [~,preds]= max(nn.a{end});
    preds = preds';

end