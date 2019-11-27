buildingDir = fullfile(toolboxdir('vision'), 'visiondata', 'building');
buildingScene = imageDatastore(buildingDir);

I1 = readimage(buildingScene, 1);
I2 = readimage(buildingScene, 2);

I1_gray = rgb2gray(I1);
I2_gray = rgb2gray(I2);

% get points
points1 = detectHarrisFeatures(I1_gray);
points2 = detectHarrisFeatures(I2_gray);

% get features
[features1, points1] = extractFeatures(I1_gray, points1);
[features2, points2] = extractFeatures(I2_gray, points2);

loc1 = points1.Location;
loc2 = points2.Location;

[match,match_fwd,match_bkwd] = match_features(double(features1.Features),double(features2.Features));

H = ransac_homography(loc1(match(:,1),:),loc2(match(:,2),:));

I = stitch(I1,I2,H);

figure()
imshow(I)

%%
function best_H = ransac_homography(p1,p2)
    thresh = sqrt(2); % threshold for inlier points
    p = 1-1e-4; % probability of RANSAC success
    w = 0.5; % fraction inliers
	
    % n: number of correspondences required to build the model (homography)
    n = 4; %need at least 4 points to construct homograhpy
    % number of iterations required
    % from the lecture given the probability of RANSAC success, and fraction of inliers
    k = log(1-p)/log(1-w^n); % k =log(1-p)/log(1-w^n))
	
    num_pts = size(p1,1);
    best_inliers = 4;
    best_H = eye(3);
    
    [np1,~] = size(p1);
    [np2,~] = size(p2);
    for iter = 1:k
        % randomly select n correspondences from p1 and p2
        % use these points to compute the homography
        random_index = randperm(np1);
        p1_sample = p1(random_index(1:n),:);
        p2_sample = p2(random_index(1:n),:);
        H = compute_homography(p1_sample, p2_sample);
	
        % transform p2 to homogeneous coordinates
        p2_h = [p2';ones(1,size(p2,1))];
        % estimate the location of correspondences given the homography
        p1_hat = H*p2_h;
        % convert to image coordinates by dividing x and y by the third coordinate
        p1_hat = p1_hat./p1_hat(3,:);
        % compute the distance between the estimated correspondence location and the 
        % putative correspondence location
        dist = pdist2(p1_hat(1:2,:).',p1);
        % inlying points have a distance less than the threshold thresh defined previously
        num_inliers = sum(sum(dist <thresh));
		
        if num_inliers > best_inliers
            best_inliers = num_inliers;
            best_H = H;
        end
    end
end

%%
function [match,match_fwd,match_bkwd] = match_features(f1,f2)
    %% INPUT
    %% f1,f2: [ number of points x number of features ]
    %% OUTPUT
    %% match, match_fwd, match_bkwd: [ indices in f1, corresponding indices in f2 ]
    % get matches using pdist2 and the ratio test with threshold of 0.7
    % fwd matching
    % calculate distance between f1 and f2, rows of f1 = feature1,2,...n;
    % columns of f1 are 64 descriptors (gradient, deviation, mean,..)
    % f1 is 758x64 matrix, f2 = 742x64 matrix -> D is 758x742 matrix where
    % D(1,1) = distance(f1(feature1),f2(feature2))
    D = pdist2(f1,f2);
    % sort those distance in ascending order, for each feature in f1, sort
    % from small to big distance from f1 to each feature of f2
    % Dsorted(1,:) = distance(f1(feature1),f2(feature(1:742)) sorted in
    % ascending order
    % Dindex = features of f2
    [Dsorted, Dindex] = sort(D,2); % sort in rows
    % for each feature in f1 find 2 features in f2 have smallest distance 
    % from it, store those distances  
    nearest_neighbor = Dsorted(:,1);
    nearest_neighbor2 = Dsorted(:,2);
    confidences_ratio = nearest_neighbor./ nearest_neighbor2;
    % find the ratio that < 0.7
    % return index of confidences_ratio (feature i of f1 that have ratio <0.7)
    accepted_feature_f1 = find(confidences_ratio < 0.7); 
    % number of accepted points:
    nb_accepted_feature_f1 = size(accepted_feature_f1,1); % nb of rows of accepted points
    match_fwd = zeros(nb_accepted_feature_f1,2);
    match_fwd(:,1) = accepted_feature_f1; % accepted features of f1
    match_fwd(:,2) = Dindex(accepted_feature_f1); % feature of f2 for accepted feature of f1
    
    % bkwd matching
    % calculate distance between f1 and f2, rows of f1 = feature1,2,...n; columns of f1 are 
    % 64 descriptors (gradient, deviation, mean,..)
    D = pdist2(f2,f1);
    % sort those distance in ascending order, for each feature in f1, sort
    % from small to big distance from f1 to each feature of f2
    % Dsorted(1,:) = distance(f1(feature1),f2(feature(1:742)) sorted in
    % ascending order
    % Dindex = features of f2
    [Dsorted, Dindex] = sort(D,2); % sort in rows
    % for each feature in f1 find 2 features in f2 have smallest distance 
    % from it, store those distances  
    nearest_neighbor = Dsorted(:,1);
    nearest_neighbor2 = Dsorted(:,2);
    confidences_ratio = nearest_neighbor./ nearest_neighbor2;
    % find the ratio that < 0.7
    % return index of confidences_ratio (feature i of f1 that have ratio <0.7)
    accepted_feature_f2 = find(confidences_ratio < 0.7);
    % number of accepted points:
    nb_accepted_feature_f2 = size(accepted_feature_f2,1); % nb of rows of accepted points
    match_bkwd = zeros(nb_accepted_feature_f2,2);
    match_bkwd(:,2) = accepted_feature_f2; % accepted features of f2
    match_bkwd(:,1) = Dindex(accepted_feature_f2); % feature of f1 for accepted feature of f2
    % fwd bkwd consistency check
    [C,D] = sortrows([match_fwd;match_bkwd]);
    E = all(C(1:end-1,:)==C(2:end,:),2);
    F = D(E);
    match = match_fwd(F,:);
end

%%
function H = compute_homography(p1,p2)		
    % use SVD to solve for H as was done in the lecture
    % Note: p1 is a 4x2 matrix, the first column of which gives the x coordinates.
    % construct matrix A: A*H = 0
    A = [p2(1,1) p2(1,2) 1 0 0 0 -p2(1,1)*p1(1,1) -p2(1,2)*p1(1,1) -p1(1,1);
        0 0 0 p2(1,1) p2(1,2) 1 -p2(1,1)*p1(1,2) -p2(1,2)*p1(1,2) -p1(1,2);
        p2(2,1) p2(2,2) 1 0 0 0 -p2(2,1)*p1(2,1) -p2(2,2)*p1(2,1) -p1(2,1);
        0 0 0 p2(2,1) p2(2,2) 1 -p2(2,1)*p1(2,2) -p2(2,2)*p1(2,2) -p1(2,2);
        p2(3,1) p2(3,2) 1 0 0 0 -p2(3,1)*p1(3,1) -p2(3,2)*p1(3,1) -p1(3,1);
        0 0 0 p2(3,1) p2(3,2) 1 -p2(3,1)*p1(3,2) -p2(3,2)*p1(3,2) -p1(3,2);
        p2(4,1) p2(4,2) 1 0 0 0 -p2(4,1)*p1(4,1) -p2(4,2)*p1(4,1) -p1(4,1);
        0 0 0 p2(4,1) p2(4,2) 1 -p2(4,1)*p1(4,2) -p2(4,2)*p1(4,2) -p1(4,2)];
    % construct U,D,V matrices from A using singular value decomposition: 
    % A = U*D*V.'; H is the last column of V;
    [U, D, V] = svd(A);
    % X = [h11 h12 ... h33].'; H = [h11 h12 h13; h21 h22 h23; h31 h32 h33];
    % X = last column of V, make h33 = 1 for normalization
    X = V(:,end)/V(end,end);
    % contruct H
    H = reshape(X,3,3).';  
end

%%
function I_ = stitch(I1,I2,H)
	[sy2,sx2,sz2] = size(I2);
	[x2,y2,z2] = meshgrid(1:sx2,1:sy2,1:sz2);
	
	% map I2 to I1
	p1_hat = H*[x2(:),y2(:),ones(numel(x2),1)]';
	p1_hat = p1_hat(1:2,:)./(repmat(p1_hat(3,:),2,1)+eps);
	
	% create new dimensions to accomodate points from I2
	p1_hat_xmax = max(p1_hat(1,:));
	p1_hat_xmin = min(p1_hat(1,:));
	p1_hat_ymax = max(p1_hat(2,:));
	p1_hat_ymin = min(p1_hat(2,:));
	
	xmin = round(floor(min(p1_hat_xmin,0)));
	xmax = round(ceil(max(p1_hat_xmax,sx2)));
	ymin = round(floor(min(p1_hat_ymin,0)));
	ymax = round(ceil(max(p1_hat_ymax,sy2)));
	
	% create images for mapping
	I1_ = uint8(zeros(ymax - ymin, xmax - xmin,3));
	I2_ = uint8(zeros(ymax - ymin, xmax - xmin,3));
	I_ = uint8(zeros(ymax - ymin, xmax - xmin,3));
	
	% I1 is just translated in I_
	I1_(1-ymin:sy2-ymin, 1-xmin:sx2-xmin,:) = I1;
	
	% map I_ to I2 (translation then homography)
	[sy2_,sx2_,sz2_] = size(I2_);
	[x2_,y2_,z2_] = meshgrid(1:sx2_,1:sy2_,1:sz2_);
	p2_hat = H\[x2_(:)+xmin, y2_(:)+ymin, ones(numel(x2_),1)]';
	p2_hat = round(p2_hat(1:2,:)./(repmat(p2_hat(3,:),2,1)+eps));
	
	% keep only the valid coordinates
	good_x = p2_hat(1,:) > 0 & p2_hat(1,:) <= sx2;
	good_y = p2_hat(2,:) > 0 & p2_hat(2,:) <= sy2;
	good = good_x & good_y; good = good(:);
	
	% valid in I_
	ind2 = sub2ind(size(I2_),y2_(good),x2_(good),z2_(good));
	% valid in I2
	ind2_hat = sub2ind(size(I2),p2_hat(2,good)',p2_hat(1,good)',z2_(good));
	
	% I2 transformed by homography in I_
	I2_(ind2) = I2(ind2_hat);
	
	% nonoverlapping regions do not require blending
	I2_sum = sum(I2_,3); I1_sum = sum(I1_,3);
	[no_blend_y,no_blend_x] = find(I2_sum == 0 | I1_sum == 0);
	
	no_blend_ind_r = sub2ind(size(I_),no_blend_y,no_blend_x,ones(size(no_blend_x)));
	no_blend_ind_g = sub2ind(size(I_),no_blend_y,no_blend_x,2*ones(size(no_blend_x)));
	no_blend_ind_b = sub2ind(size(I_),no_blend_y,no_blend_x,3*ones(size(no_blend_x)));
	
	I_(no_blend_ind_r) = I2_(no_blend_ind_r) + I1_(no_blend_ind_r);
	I_(no_blend_ind_g) = I2_(no_blend_ind_g) + I1_(no_blend_ind_g);
	I_(no_blend_ind_b) = I2_(no_blend_ind_b) + I1_(no_blend_ind_b);
	
	% overlapping regions require blending		
	[blend_y,blend_x] = find(I2_sum > 0 & I1_sum > 0);
	
	blend_ind_r = sub2ind(size(I_),blend_y,blend_x,ones(size(blend_x)));
	blend_ind_g = sub2ind(size(I_),blend_y,blend_x,2*ones(size(blend_x)));
	blend_ind_b = sub2ind(size(I_),blend_y,blend_x,3*ones(size(blend_x)));
	
	I_(blend_ind_r) = uint8(double(I2_(blend_ind_r))*.5 + double(I1_(blend_ind_r))*.5);
	I_(blend_ind_g) = uint8(double(I2_(blend_ind_g))*.5 + double(I1_(blend_ind_g))*.5);
	I_(blend_ind_b) = uint8(double(I2_(blend_ind_b))*.5 + double(I1_(blend_ind_b))*.5);
end

%%
function plot_corr(I1,I2,p1,p2)
	I = [I1,I2];
	[sy,sx] = size(I1);
	
	figure()
	imshow(I)
	hold on;
	plot(p1(:,1),p1(:,2),'bo');
	plot(sx + p2(:,1),p2(:,2),'rx');
	plot([p1(:,1),sx+p2(:,1)]',[p1(:,2),p2(:,2)]','g-');
	hold off;
end