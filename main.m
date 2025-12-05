clear; clc; 
close all;

% % Parameters of RANSAC via fundamental matrix
parameters.minPtNum = 6;    % minimal number for model fitting
parameters.iterNum = 2000;  % maximum number of trials
parameters.thDist = 0.1;   % distance threshold for inliers

datasetname='PTS';
pathname=strcat('Imgs/', datasetname, '/');
outpath=strcat(pathname, 'results/');
imgs_format='*.jpg'; 
dir_folder=dir(strcat(pathname, imgs_format));
if isempty(dir_folder)
    error('no input images found in this folder');
end

path1=sprintf('%s%s',pathname,dir_folder(1).name); % target image
path2=sprintf('%s%s',pathname,dir_folder(2).name); % reference image

depth_img1 = double(rgb2gray(imread([pathname, dir_folder(1).name(1:end-4), '_depth_vitl.png'])));
depth_img1 = 1e1./depth_img1;
depth_img1(depth_img1==inf)=1e4;

depth_img2 = double(rgb2gray(imread([pathname, dir_folder(2).name(1:end-4), '_depth_vitl.png'])));
depth_img2 = 1e1./depth_img2;
depth_img2(depth_img2==inf)=1e4;

%% Read images.
%-------------
img1 = im2double(imread(path1));   img2 = im2double(imread(path2));    
%% detect and match sift features and line features
[pts1, pts2] = siftMatch(img1, img2); 

%% depth-based RANSAC
[matches1, matches2]=fundRANSAC(pts1, pts2, depth_img1, parameters); % feature ransac via fundamental matrix

%% epipolar geometry estimation
[H_inf1, fund_e1] = calcFundMatrixofHe(matches1, matches2, depth_img1);
[H_inf2, fund_e2] = calcFundMatrixofHe(matches2, matches1, depth_img2);

%% optimal image warping
[warped_img1, warped_img2, panorama, quality_errors]=optimalBackwardMapping(img1,img2,depth_img1,depth_img2,H_inf1,fund_e1,H_inf2, fund_e2);
