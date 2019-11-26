function K = KernelIntersection(X1, X2)
   % computes a histogram intersection kernel
   %
   % Input:
   % - X1: an n x d dimensional feature matrix where n is the number of observations, and d is the number of features.
   % - X2: an m x d dimensional feature matrix where m is the number of observations, and d is the number of features.
   % Output:
   % - K: an n x m dimensional kernel matrix where K(i,j) stores a histogram intersection between the data point i in X1, and data point j in X2.
   K = zeros(size(X1,1),size(X2,1));
   for rows = 1: size(X1,1)
       for cols = 1:size(X2,1)
           K(rows,cols) = sum(min(X1(rows,:),X2(cols,:)));
       end
   end
end