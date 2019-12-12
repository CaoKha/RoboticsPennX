N = 20;
Ic = ImageCarvingFunction(N);
function Ic = ImageCarvingFunction(N)

% N: number of vertical seams you have to remove


% read image
I = im2double(imread('waterfall.png'));
% get grayscale image
Ig0 = rgb2gray(im2double(I));

% colored image 
Ic = I; 


for iIter = 1:N

    Ig = rgb2gray(Ic);
    Gx = imfilter(Ig,.5*[-1 0 1],'replicate');
    Gy = imfilter(Ig,.5*[-1 0 1]','replicate');
    E = abs(Gx) +  abs(Gy); % energy
    
    
    M = zeros(size(Ig)); % cumulative energy
    
    % your CODE starts here
    [rows, columns] = size(E);
    P = zeros(rows,columns); % path function
    
    % set first rows of M to first rows of E 
    M(1,:) = E(1,:);
    
    % construct cumulative energy map M and path matrix P
    for i = 2: rows
        for j = 2: columns-1
            [minimum,index_min] = min([M(i-1,j-1), M(i-1,j), M(i-1,j+1)]);
            M(i,j) = E(i,j) + minimum;
            P(i,j) = index_min -2;
        end
        [minimum_col1,index_col1] = min([M(i-1,1), M(i-1,2)]);
        M(i,1) = E(i,1) + minimum_col1;
        P(i,1) = index_col1 - 1;
        [minimum_last_col,index_last_col] = min([M(i-1,columns-1), M(i-1,columns)]);
        M(i,columns) = E(i,columns) + minimum_last_col;
        P(i,columns) = index_last_col - 2;
    end
    
    % construct vertical seam 
    verticalSeam = zeros(rows,1);
    [~,col_index] =  min(M(rows,:));
    verticalSeam(rows) = col_index;
    
    for i = (rows-1):-1:1
        verticalSeam(i) = verticalSeam(i+1) + P(i+1,verticalSeam(i+1));
    end
    
    % reduce width of image by removing the vertical seam
    reduceIm = zeros(rows,columns-1, 3);
    
    for i = 1:rows
        reduceIm(i, 1:(verticalSeam(i)-1),:) = Ic(i,1:(verticalSeam(i)-1),:);
        reduceIm(i,verticalSeam(i):(columns-1),:) = Ic(i,(verticalSeam(i)+1):columns,:);
    end
    
    Ic = reduceIm;
    % your CODE ends here


end
figure(1),imshow(I);
figure(2),imshow(Ic);
figure(3),imshow(E);


end