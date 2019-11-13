warp_frac = 1/2;
dissolve_frac = 1/2;
M = ImageMorphingTriangulation(warp_frac,dissolve_frac);

function M = ImageMorphingTriangulation(warp_frac,dissolve_frac)

if nargin < 1
    warp_frac = .5;
end

if nargin < 2
    dissolve_frac= warp_frac; 
end


% ream images
I = im2double(imread('a.png'));
J = im2double(imread('c.png'));

% load mat file with points, variables Ip,Jp
 load('points.mat');
 
% initialize ouput image (morphed)
M = zeros(size(I));

%  Triangulation (on the mean shape)
MeanShape = (1/2)*Ip+(1/2)*Jp;
TRI = delaunay(MeanShape(:,1),MeanShape(:,2));


% number of triangles
TriangleNum = size(TRI,1); 

% find coordinates in images I and J
CordInI = zeros(3,3,TriangleNum);
CordInJ = zeros(3,3,TriangleNum);

for i =1:TriangleNum
  for j=1:3
    
    CordInI(:,j,i) = [ Ip(TRI(i,j),:)'; 1];
    CordInJ(:,j,i) = [ Jp(TRI(i,j),:)'; 1]; 
    
  end
end

% create new intermediate shape according to warp_frac
Mp = (1-warp_frac)*Ip+warp_frac*Jp; 

 
% create a grid for the morphed image
[x,y] = meshgrid(1:size(M,2),1:size(M,1));

% for each element of the grid of the morphed image, find  which triangle it falls in
TM = tsearchn([Mp(:,1) Mp(:,2)],TRI,[x(:) y(:)]);


% YOUR CODE STARTS HERE
% suppress NaN components in TM, x and y
idx = ~isnan(TM); x = x(idx); y = y(idx); TM = TM(idx);

% Calculate matrix 3 corners triangle of intermedate image
CordInM = zeros(3,3,TriangleNum);
for i =1:TriangleNum
  for j=1:3  
    CordInM(:,j,i) = [ Mp(TRI(i,j),:)'; 1]; 
  end
end

% compute barycentric coordinates of intermedate image
barycentric = zeros(size(x,1),3);
PointI = zeros(size(x,1),3);
PointJ = zeros(size(x,1),3);
PointM = zeros(size(x,1),3);
for k = 1:size(x)
    barycentric(k,1:3) = CordInM(:,:,TM(k))\[x(k);y(k);1];
    PointI(k,1:3) = (CordInI(:,:,TM(k))*barycentric(k,1:3)')';
    PointJ(k,1:3) = (CordInJ(:,:,TM(k))*barycentric(k,1:3)')';
    PointM(k,1:3) = (CordInM(:,:,TM(k))*barycentric(k,1:3)')';
end

IndI = round(PointI(:,1:2)./PointI(:,3));
IndJ = round(PointJ(:,1:2)./PointJ(:,3));
IndM = round(PointM(:,1:2)./PointM(:,3));


IndI = sub2ind(size(I),IndI(:,2),IndI(:,1));
IndJ = sub2ind(size(I),IndJ(:,2),IndJ(:,1));
IndM = sub2ind(size(I),IndM(:,2),IndM(:,1));
% YOUR CODE ENDS HERE



% cross-dissolve
M(IndM)=(1-dissolve_frac)* I(IndI)+ dissolve_frac * J(IndJ);
M(IndM+53200)=(1-dissolve_frac)* I(IndI+53200)+ dissolve_frac * J(IndJ+53200);
M(IndM+2*53200)=(1-dissolve_frac)* I(IndI+2*53200)+ dissolve_frac * J(IndJ+2*53200);

figure(100);
subplot(1,3,1);
imshow(I);
hold on;
triplot(TRI,Ip(:,1),Ip(:,2))
hold off;
title('First')

subplot(1,3,2);
imshow(M);
hold on;
triplot(TRI,Jp(:,1),Jp(:,2))
hold off
title('Morphed')

subplot(1,3,3);
imshow(J);
hold on;
triplot(TRI,Jp(:,1),Jp(:,2))
hold off
title('Second')

end