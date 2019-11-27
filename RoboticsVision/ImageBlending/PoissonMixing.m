I = PoissonMixingGradients;
function I = PoissonMixingGradients
    
% read images 
target= im2double(imread('target_2.jpg')); 
source= im2double(imread('source_2.jpg')); 
mask=imread('mask_2.bmp');

row_offset=130;
col_offset=10; 

source_scale=0.6;

source =imresize(source,source_scale);
mask =imresize(mask,source_scale);

size_target=size(target);
size_source=size(source);

N=sum(mask(:)); % N: Number of unknown pixels == variables

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
    
    Np(ib)=  double((row_offset+i> 1))+ ...
             double((col_offset+j> 1))+ ...
             double((row_offset+i< size(target,1))) + ...
             double((col_offset+j< size(target,2)));
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
          if abs(source(i,j,color)-source(i-1,j,color)) >= abs((target(row_offset+i,col_offset+j,color)-target(row_offset+i-1,col_offset+j,color)))
              b(ib)=b(ib)+ target(row_offset+i-1,col_offset+j,color)*(1-mask(i-1,j))+...
                          source(i,j,color)-source(i-1,j,color);
          else
              b(ib)=b(ib)+ target(row_offset+i-1,col_offset+j,color)*(1-mask(i-1,j))+...
                          (target(row_offset+i,col_offset+j,color)-target(row_offset+i-1,col_offset+j,color));
          end
      end

      if (i<size(mask,1))
          if abs(source(i,j,color)-source(i+1,j,color)) >= abs((target(row_offset+i,col_offset+j,color)-target(row_offset+i+1,col_offset+j,color)))
              b(ib)=b(ib)+  target(row_offset+i+1,col_offset+j,color)*(1-mask(i+1,j))+ ...
                           source(i,j,color)-source(i+1,j,color);
          else
              b(ib)=b(ib)+  target(row_offset+i+1,col_offset+j,color)*(1-mask(i+1,j))+ ...
                           (target(row_offset+i,col_offset+j,color)-target(row_offset+i+1,col_offset+j,color));
          end      
      end

      if (j>1)
          if abs((source(i,j,color)-source(i,j-1,color))) >= abs((target(row_offset+i,col_offset+j,color)-target(row_offset+i,col_offset+j-1,color)))
              b(ib)= b(ib) +  target(row_offset+i,col_offset+j-1,color)*(1-mask(i,j-1))+...
                           source(i,j,color)-source(i,j-1,color);
          else
              b(ib)= b(ib) +  target(row_offset+i,col_offset+j-1,color)*(1-mask(i,j-1))+...
                           (target(row_offset+i,col_offset+j,color)-target(row_offset+i,col_offset+j-1,color));
          end
      end

      if (j<size(mask,2))
          if abs((source(i,j,color)-source(i,j+1,color))) >= abs((target(row_offset+i,col_offset+j,color)-target(row_offset+i,col_offset+j+1,color)))
              b(ib)= b(ib)+ target(row_offset+i,col_offset+j+1,color)*(1-mask(i,j+1))+...
                         source(i,j,color)-source(i,j+1,color); 
          else
              b(ib)= b(ib)+ target(row_offset+i,col_offset+j+1,color)*(1-mask(i,j+1))+...
                         (target(row_offset+i,col_offset+j,color)-target(row_offset+i,col_offset+j+1,color)); 
          end
      end     

    end

    
    % solve linear system A*x = b;
    x = A\b;

    % impaint target image
     for ib=1:N
           I(row_offset+ir(ib),col_offset+ic(ib),color) = x(ib);
     end
    
end
% YOUR CODE ENDS HERE
figure(1); imshow(source);
figure(2); imshow(target);
figure(3); imshow(I);
figure(4); imshow(I - target);


end