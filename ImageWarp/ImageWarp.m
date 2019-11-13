v = VideoReader('shuttle.avi');
fr1 = readFrame(v);
im1t = im2double(rgb2gray(fr1));

hasFrame(v);
fr2 = readFrame(v); fr2 = readFrame(v);
im2t = im2double(rgb2gray(fr2));

[u,v] = multiscale_flow(im1t,im2t);

figure()
subplot(221)
imshow(im1t)
subplot(222)
imshow(im2t)
subplot(223)
imagesc(u)
subplot(224)
imagesc(v)

function [u,v] = estimate_flow(I1,I2,wsize)
    %% INPUT:
    %% I1, I2: nxm sequential frames of a video
    %% wsize: (wsize*2 + 1)^2 is the size of the neighborhood used for displacement estimation
    %% OUTPUT:
    %% u,v: nxm flow estimates in the x and y directions respectively
    
    % Compute the image gradients for the second image
    [I2_x,I2_y] = grad2d(I2);
    % The temporal gradient is the smoothed difference image
    I2_t = gauss_blur(I1-I2);
    % initialize x and y displacement to zero
    u = zeros(size(I2));
    v = zeros(size(I2));

    % loop over all pixels in the allowable range
    for i = wsize+1:size(I2_x,1)-wsize
       for j = wsize+1:size(I2_x,2)-wsize

          % Select the appropriate window
          Ix = I2_x(i-wsize:i+wsize,j-wsize:j+wsize);
          Iy = I2_y(i-wsize:i+wsize,j-wsize:j+wsize);
          It = I2_t(i-wsize:i+wsize,j-wsize:j+wsize);
   
   	  d = estimate_displacement(Ix,Iy,It);

          u(i,j) = d(1,:);
          v(i,j) = d(2,:);
       end
    end
    % use medifilt2 with a 5x5 filter to reduce outliers in the flow estimate
    u = medfilt2(u,[5 5]);
    v = medfilt2(v,[5 5]);
end

function d = estimate_displacement(Ix,Iy,It)
    %% INPUT:
    %% Ix, Iy, It: m x m matrices, gradient in the x, y and t directions
    %% Note: gradient in the t direction is the image difference
    %% OUTPUT:
    %% d: least squares solution
    
    b = [ Ix(:) Iy(:) ]' * It(:);
    % to help mitigate effects of degenerate solutions add eye(2)*eps to the 2x2 matrix A
    A = [ Ix(:) Iy(:) ]' * [ Ix(:) Iy(:) ] + eye(2)*eps; 
    % use pinv(A)*b to compute the least squares solution
    d = pinv(A)*b;
end

function [I_x,I_y] = grad2d(img)
	%% compute image gradients in the x direction
	%% convolve the image with the derivative filter from the lecture
	%% using the conv2 function and the 'same' option
	dx_filter = [1/2 0 -1/2];
    I_x = conv2(img,dx_filter,'same'); 

	%% compute image gradients in the y direction
	%% convolve the image with the derivative filter from the lecture
	%% using the conv2 function and the 'same' option
	dy_filter = dx_filter.';
    I_y = conv2(img,dy_filter,'same'); 
end

function smooth = gauss_blur(img)
    %% Since the Gaussian filter is separable in x and y we can perform Gaussian smoothing by
    %% convolving the input image with a 1D Gaussian filter in the x direction then  
    %% convolving the output of this operation with the same 1D Gaussian filter in the y direction.

    %% Gaussian filter of size 5
    %% the Gaussian function is defined f(x) = 1/(sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma))
    x = -2:2;
    sigma = 1;
    gauss_filter = 1/(sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma));

    %% using the conv2 function and the 'same' option
    %% convolve the input image with the Gaussian filter in the x
    smooth_x = conv2(img,gauss_filter,'same');
    %% convolve smooth_x with the transpose of the Gaussian filter
    smooth = conv2(smooth_x,gauss_filter.','same');
end

function [u,v] = multiscale_down(I1,I2,u,v,l)
    % If the base pyramid level has been reached, return.
    if l == 0
        return
    end
    % Otherwise, upsample the previous flow estimate by a factor of 2 using imresize with 
    % bicubic interpolation. The flow values should be doubled.
    u = 2*imresize(u,2,'bicubic');
    v = 2*imresize(v,2,'bicubic');
    % Warp the input image, I2, according to the upsampled flow estimate.
    I2_ = warp_image(I2,u,v);
    % Estimate the incremental flow update using your estimate_flow function and the warped 
    % input image with 2 as the wsize parameter.
    [u_,v_] = estimate_flow(I1,I2_,2);
    % Update the flow estimate by adding the incremental estimate above to the previous estimate.
    u = u + u_;
    v = v + v_;
end

function [u,v] = multiscale_aux(I1,I2,l,lmax)
    % Downsample the images by half using imresize and bicubic interpolation.
    % Use your gauss_blur function to smooth the result.
    I1_ = gauss_blur(imresize(I1,0.5,'bicubic'));
    I2_ = gauss_blur(imresize(I2,0.5,'bicubic'));
    % If the highest pyramid level has been reached, estimate the optical flow
    % on the downsampled images with your estimate_flow function using 2 as the wsize parameter.
    if l == lmax
    	[u,v] = estimate_flow(I1_,I2_,2);
    % If we are beyond the highest pyramid level, estimate the optical flow
    % on the input images (not the downsampled images) with your estimate_flow function 
    % using 2 as the wsize parameter.
    elseif l > lmax
        l = lmax;
        [u,v] = estimate_flow(I1,I2,2);
    % Otherwise, increment the current level and continue up the pyramid (i.e. recurse)
    % using the downsampled images.
    else
        [u,v] = multiscale_aux(I1_,I2_,l+1,lmax);
    end
    % After flow has been estimated at the current level, pass this estimate along with
    %   (not the downsampled images) to multiscale_down for iterative 
    % improvement of the flow estimate.
    [u,v] = multiscale_down(I1,I2,u,v,l);
end

function [u,v] = multiscale_flow(I1,I2)
    % The number of pyramid levels will be determined by the image size
    % At the highest pyramid level the smallest image dimension will be around 
    % 30 pixels.
    lmax = round(log2(min(size(I1))/30));
    % The pyramidal approach can be implemented with a recursive strategy
    [u,v] = multiscale_aux(I1,I2,1,lmax);
end

function warp = warp_image(I,u,v)
    %% INPUT:
    %% I: image to be warped
    %% u,v: x and y displacement
    %% OUTPUT:
    %% warp: image I deformed by u,v
    
    % initialize warp as zeros
    warp = zeros(size(I));
    % construct warp so that warp(x,y) = I(x + u, y + v)
    [X,Y] = meshgrid(1:size(I,2),1:size(I,1));
    warpedX = X + round(u);
    warpedY = Y + round(v);
    % interpolate each point by using "interp2" function
    warp = interp2(X,Y,I,warpedX,warpedY,'cubic',0);
end

function [I1, I2, I3, I4] = test_images()
	I1 = repmat(1:20,20,1);
	I2 = [0, 0, 0, 1:17 ];
	I2 = repmat(I2,20,1);

	I3 = I1 + I1';
	I4 = I2 + I2';
end