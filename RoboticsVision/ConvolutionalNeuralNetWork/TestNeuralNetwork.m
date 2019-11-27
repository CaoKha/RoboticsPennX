load('X.mat');
load('Y.mat');   
X=StandardizeData(X);
[train_x, train_y, test_x, test_y] = splitData(X, Y);
nn=TrainNetwork(train_x,train_y);
preds=PredictNetwork(nn,test_x);

% check prediction precision
ground_truth = zeros(size(preds,1),1);
for i = 1:size(preds,1)
    for j = 1:size(test_y,2)
        if test_y(i,j) == 1
            ground_truth(i) = j;
        end
    end
end
count = 0;
for id = 1: size(preds,1)
    if preds(id) == ground_truth(id)
        count = count +1;
    end
end
prediction_precision = count/size(preds,1);