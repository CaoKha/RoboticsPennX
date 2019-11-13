clear all;
close all;
img = imread('peppers.png');
img_gray = double(rgb2gray(img));

img_gray_smooth = gauss_blur(img_gray);
[I_x,I_y] = grad2d(img_gray_smooth);

I_xx = gauss_blur(I_x.*I_x);
I_yy = gauss_blur(I_y.*I_y);
I_xy = gauss_blur(I_x.*I_y);

k = 0.06;
%% Use the corner score equation from the lecture.
R = (I_xx.*I_yy-I_xy.^2) - k*(I_xx + I_yy).^2;

%%
r = 5;
thresh = 10000;
hc = nmsup(R,r,thresh);

figure()
imshow(img)
hold on;
plot(hc(:,1), hc(:,2), 'rx')
hold off;
%%
function smooth = gauss_blur(img)
%% Since the Gaussian filter is separable in x and y we can perform Gaussian smoothing by
%% convolving the input image with a 1D Gaussian filter in the x direction then
%% convolving the output of this operation with the same 1D Gaussian filter in the y direction.

%% Gaussian filter of size 5
%% the Gaussian function is defined f(x) = 1/(sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2))
x = -2:2;
sigma = 1;
gauss_filter = 1/(sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2));

%% using the conv2 function and the 'same' option
%% convolve the input image with the Gaussian filter in the x
smooth_x = conv2(img,gauss_filter,'same');
%% convolve smooth_x with the transpose of the Gaussian filter
smooth = conv2(smooth_x,gauss_filter.','same');
end

%%
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

function loc = nmsup(R,r,thresh)
%% Step 1-2 must be performed in a way that allows you to
%% preserve location information for each corner.
[sy,sx] = size(R);
[x,y] = meshgrid(1:sx,1:sy);

%% Step 1: eliminate values below the specified threshold.
Rcorner = zeros(1);
R_Xindex = zeros(1);
R_Yindex = zeros(1);
count = 1;
for i = 1:sy
    for j = 1:sx
        if R(i,j) < thresh
            R(i,j) = 0;
        else
            Rcorner(count) = R(i,j);
            R_Xindex(count) = j;
            R_Yindex(count) = i;
            count = count + 1;
        end
    end
end
%% Step 2: Sort the remaining values in decreasing order.
[Rcorner, Rcorner_index] = sort(Rcorner,'descend');

%% Step 3: Starting with the highest scoring corner value, if
%% there are corners within its r neighborhood remove
%% them since their scores are lower than that of the corner currently
%% considered. This is true since the corners are sorted
%% according to their score and in decreasing order.
max_element = size(Rcorner,2);
for k = 1: max_element
    for u = -r:r
        for v = -r:r
            if (u^2 + v^2 <= r^2) && (R_Yindex(Rcorner_index(k)) + u <= sy)...
                    && (R_Yindex(Rcorner_index(k)) + u > 0)...
                    && (R_Xindex(Rcorner_index(k)) + v > 0)...
                    && (R_Xindex(Rcorner_index(k)) + v <= sx)...
                    && (R(R_Yindex(Rcorner_index(k)) + u,R_Xindex(Rcorner_index(k)) + v) <=...
                    R(R_Yindex(Rcorner_index(k)),R_Xindex(Rcorner_index(k))))... 
                    && (u ~= 0 || v ~= 0)
                R(R_Yindex(Rcorner_index(k)) + u,R_Xindex(Rcorner_index(k)) + v) = 0;
            end
        end
    end
end

Rcorner2 = zeros(1);
R_Xindex2 = zeros(1);
R_Yindex2 = zeros(1);
count2 = 1;
for i = 1:sy
    for j = 1:sx
        if R(i,j) < thresh
            R(i,j) = 0;
        else
            Rcorner2(count2) = R(i,j);
            R_Xindex2(count2) = j;
            R_Yindex2(count2) = i;
            count2 = count2 + 1;
        end
    end
end

[Rcorner2, Rcorner_index2] = sort(Rcorner2,'descend');
%% The variable loc should contain the sorted corner locations which
%% survive thresholding and non-maximum suppression with
%% size(loc): nx2
%% loc(:,1): x location
%% loc(:,2): y location
max_element2 = size(Rcorner2,2);
loc = zeros(max_element2,2);
for k = 1: max_element2
    loc(k,1) = R_Xindex2(Rcorner_index2(k));
    loc(k,2) = R_Yindex2(Rcorner_index2(k));
end

end
