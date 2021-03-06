function D = derivative_matrix(n)
    syms k
    % your implementation...
    D = sym(zeros(n,8));
    for i = 1:n
        for j = 1:8
            if i == 1
                D(i,j) = k^(j-1);
            else
                D(i,j) = diff(D(i-1,j));
            end
        end
    end
            
end