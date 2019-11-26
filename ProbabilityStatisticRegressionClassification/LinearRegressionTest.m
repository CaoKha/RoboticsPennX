load('CarData.mat');
w=LinearRegression(trainsetX,trainsetY);
y_hat = Prediction(testsetX,w);
err=Error(y_hat,testsetY);

function w = LinearRegression(X,y)
    % learn a linear regression model
    %
    % Input:
    % - X: n x 7 dimensional feature matrix where every row depicts a particular observation and every column denotes a feature.
    % - y: n x 1 dimensional ground truth MPG labels
    % Output:
    % - w: 7 x 1 dimensional weights of a learned linear regression model
    
    w= (X'*X)\(X'*y);
end

function y_hat = Prediction(X,w)
    % make predictions using a learned linear regression model
    %
    % Input:
    % - X: n x 7 dimensional feature matrix where every row depicts a particular observation and every column denotes a feature.
    % - w: 7 x 1 dimensional weights of a learned linear regression model 
    % Output:
    % - y_hat: n x 1 dimensional MPG predictions of our model
    
    y_hat= X*w;
end

function err = Error(y_hat,y)
    % estimate the error of our model
    %
    % Input:
    % - y_hat: n x 1 dimensional MPG predictions of our model
    % - y: n x 1 dimensional ground truth MPG values
    % Output:
    % - err: average L2 error of our predictions with respect to the ground truth values
    err= 1/(size(y,1))*1/2*sum((y -y_hat).^2);
end