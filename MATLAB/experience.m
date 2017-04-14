clear;
close all;
%% image data convert to Graph Matrix

% params
RESIZE = 50;
r = 5;
sigma_i = 0.1;
sigma_x = 4.0;
addpath('./data/input/')
[graph resize_img] = img2graph('37073.jpg', RESIZE, r, sigma_i, sigma_x);

% Graph Matrix
W = graph;
D = diag(sum(W,2));
L = D-W;

%% Normalized Cut
k = 2;
ncut_clusters = ncut(L,D,k);
image_seg_ncut = reshape(ncut_clusters,RESIZE,RESIZE);
image_seg_ncut(image_seg_ncut ==1) = 0;
image_seg_ncut(image_seg_ncut ==2) = 100;

%% Fast Normalized Cut with linear constraints
A = D^(-0.5)*W*D^(-0.5);
v = fast_ncut(A, 1, 10000, 10^(-9));

image_seg_fast_ncut = kmeans(v,2);
image_seg_fast_ncut = reshape(image_seg_fast_ncut,RESIZE,RESIZE);
% Segmentation result monocrize
image_seg_fast_ncut(image_seg_fast_ncut == 2) = 100;
image_seg_fast_ncut(image_seg_fast_ncut == 1) = 0;

%% Show result
figure;
subplot(1,3,1);
image(resize_img);
title('Original image')
subplot(1,3,2);
image(image_seg_ncut);
title('Normalized Cut')
subplot(1,3,3);
image(image_seg_fast_ncut);
title('Fast Ncut')

colormap(gray)
