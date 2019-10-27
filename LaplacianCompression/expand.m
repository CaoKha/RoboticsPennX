function g = expand(I)

    % Input:
    % I: the input image
    % Output:
    % g: the image after the expand operation

    % Please follow the instructions to fill in the missing commands.
    
    % 1) Create the expanded image. 
    % The new image should be twice the size of the original image.
    % So, for an n x n image you will create an empty 2n x 2n image
    % Fill every second row and column with the rows and columns of the original image
    % i.e., 1st row of I -> 1st row of expanded image
    %       2nd row of I -> 3rd row of expanded image
    %       3rd row of I -> 5th row of expanded image, and so on
    [rows,cols,dim] = size(I);
    Iexpanded = zeros(2*rows,2*cols,dim);
    u = 1; v = 1;
    for color = 1:dim
        for i = 1:rows
            for j = 1:cols
                Iexpanded(u,v,color)=I(i,j,color);
                v = v + 2;
            end
            v = 1;
            u = u + 2;
        end
        u = 1;
    end
    % 2) Create a Gaussian kernel of size 5x5 and 
    % standard deviation equal to 1 (MATLAB command fspecial)
    Gfil = fspecial('gaussian',[5 5],1);
    % 3) Convolve the input image with the filter kernel (MATLAB command imfilter)
    % Tip: Use the default settings of imfilter
    % Remember to multiply the output of the filtering with a factor of 4
    FilteredIm = 4*imfilter(Iexpanded,Gfil);
    g = FilteredIm;
    

end