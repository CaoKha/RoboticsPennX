load('ImageDataTrain.mat');
Xtrain=StandardizeData(data.trainX); 
Ytrain=data.trainY;

load('ImageDataTest.mat');
Xtest=StandardizeData(data.testX); 
Ytest=data.testY;

model = fitcsvm(Xtrain,Ytrain,'KernelFunction','KernelIntersection');
[preds,~] = predict(model,Xtest);

function X_new=StandardizeData(X)
   % standardizes data
   %
   % Input:
   % - X: an n x d dimensional feature matrix where n is the number of observations, and d is the number of features.
   % Output:
   % - X_new: an n x d dimensional normalized feature matrix.
   % calculate mean
   mean = zeros(size(X,2));
   for cols = 1:size(X,2)
       mean(cols) = sum(X(1:size(X,1),cols))/size(X,1);
   end
   %calculate standard deviation
   sd = zeros(size(X,2));
   for cols = 1:size(X,2)
       sd(cols) = sqrt(sum((X(1:size(X,1),cols) - mean(cols)).^2)/(size(X,1)-1));
   end
   
   X_new = zeros(size(X,1),size(X,2));
   for rows = 1:size(X,1)
       for cols = 1:size(X,2)
           X_new(rows,cols) = (X(rows,cols) - mean(cols)) / sd(cols);
       end
   end  
end

