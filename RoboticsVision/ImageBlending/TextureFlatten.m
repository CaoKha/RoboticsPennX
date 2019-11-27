I = TextureFlattening;

function I = TextureFlattening
    
% read image and mask
target = im2double(imread('bean.jpg')); 
mask = imread('mask_bean.bmp');

% edge detection
Edges = edge(rgb2gray(target),'canny',0.1);



N=sum(mask(:));  % N: Number of unknown pixels == variables



% YOUR CODE STARTS HERE
% enumerating pixels in the mask
mask_id = zeros(size(mask));
mask_id(mask) = 1:N; 

% neighborhood size for each pixel in the mask
[ir,ic] = find(mask);

Np = zeros(N,1); 

for ib=1:N
    
    i = ir(ib);
    j = ic(ib);
    
    Np(ib)=  double((i> 1))+ ...
             double((j> 1))+ ...
             double((i< size(target,1))) + ...
             double((j< size(target,2)));
end

% compute matrix A
A = sparse(N,N);
for ib = 1:N
    i = ir(ib);
    j = ic(ib);
    A(ib,ib) = Np(ib);
    if (i>1) && (mask(i-1,j) == 1) 
        A(ib,mask_id(i-1,j)) = -1;
    end
    if (i<size(mask,1)) && (mask(i+1,j) == 1) 
        A(ib,mask_id(i+1,j)) = -1;
    end
    if (j>1) && (mask(i,j-1) == 1) 
        A(ib,mask_id(i,j-1)) = -1;
    end
    if (j<size(mask,2)) && (mask(i,j+1) == 1) 
        A(ib,mask_id(i,j+1)) = -1;
    end
end

I = target; 

for color=1:3 % solve for each colorchannel

    % compute b for each color
    b=zeros(N,1);
    
    for ib=1:N
    
    i = ir(ib);
    j = ic(ib);
    
            
      if (i>1) 
          if (Edges(i,j) == 1) || (Edges(i-1,j) == 1)
              b(ib)=b(ib)+ target(i-1,j,color)*(1-mask(i-1,j))+...
                          target(i,j,color)-target(i-1,j,color);
          else
              b(ib)=b(ib)+ target(i-1,j,color)*(1-mask(i-1,j));                          
          end
      end

      if (i<size(mask,1))
          if (Edges(i,j) == 1) || (Edges(i+1,j) == 1)
              b(ib)=b(ib)+  target(i+1,j,color)*(1-mask(i+1,j))+ ...
                           target(i,j,color)-target(i+1,j,color);
          else
              b(ib)=b(ib)+  target(i+1,j,color)*(1-mask(i+1,j));
          end      
      end

      if (j>1)
          if (Edges(i,j) == 1) || (Edges(i,j-1) == 1)
              b(ib)= b(ib) +  target(i,j-1,color)*(1-mask(i,j-1))+...
                           target(i,j,color)-target(i,j-1,color);
          else
              b(ib)= b(ib) +  target(i,j-1,color)*(1-mask(i,j-1));
          end
      end

      if (j<size(mask,2))
          if (Edges(i,j) == 1) || (Edges(i,j+1) == 1)
              b(ib)= b(ib)+ target(i,j+1,color)*(1-mask(i,j+1))+...
                         target(i,j,color)-target(i,j+1,color); 
          else
              b(ib)= b(ib)+ target(i,j+1,color)*(1-mask(i,j+1));
          end
      end     

    end

    
    % solve linear system A*x = b;
    x = A\b;

    % impaint target image
     for ib=1:N
           I(ir(ib),ic(ib),color) = x(ib);
     end
    
end
% YOUR CODE ENDS HERE
figure(1); imshow(Edges);
figure(2); imshow(target);
figure(3); imshow(I);
figure(4); imshow(I - target);
end