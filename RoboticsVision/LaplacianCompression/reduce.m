function g = reduce(I)

    % Input:
    % I: the input image
    % Output:
    % g: the image after Gaussian blurring and subsampling

    % Please follow the instructions to fill in the missing commands.
    
    % 1) Create a Gaussian kernel of size 5x5 and 
    % standard deviation equal to 1 (MATLAB command fspecial)
    Gfil = fspecial('gaussian',[5 5],1);
    % 2) Convolve the input image with the filter kernel (MATLAB command imfilter)
    % Tip: Use the default settings of imfilter
    FilteredIm = imfilter(I,Gfil);
    % 3) Subsample the image by a factor of 2
    % i.e., keep only 1st, 3rd, 5th, .. rows and columns
    [rows,cols,dim] = size(FilteredIm);
    Isubsampled = zeros(1,1,dim);
    u = 1; v = 1;
    for color = 1:dim
        for i = 1:2:rows
            for j = 1:2:cols
                Isubsampled(u,v,color)=FilteredIm(i,j,color);
                v = v + 1;
            end
            v = 1;
            u = u + 1;
        end
        u = 1;
    end
    g = Isubsampled;
end
