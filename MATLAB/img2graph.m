function [graph, resize_img] = img2graph(file_name, RESIZE, r, sigma_i, sigma_x)
%% img2graph Image data convert to Undirected Weighted Graph.
% Inputs
% - file_name: Image data file name
% - RESIZE: resize param
% - r : threshold value if computing each node
% - sigma_i : variance of feature in Gaussian Kernel
% - sigma_x : variance of Distance in Gaussian Kernel
% Outputs
% - graph : Graph Matrix
% - resize_image : resized image data

fprintf('Image data Convert to Graph ...  \n')

%% img2graph
if nargin == 0
    file_name = '189080.jpg';
    RESIZE = 40;
    r = 10;
    sigma_i = 0.1;
    sigma_x = 4.0;
end

img = imread(file_name);
img = imresize(img,[RESIZE RESIZE]);
resize_img = img;

% Image matrix to One column data.(each column combine one column)
img =  double(img(:));
DIM = RESIZE * RESIZE;
graph = zeros(DIM,DIM);

for i = 1:DIM
    x_i = [floor((i-1)/RESIZE) mod(i-1,RESIZE)];
    for j = 1:DIM
        x_j = [floor((i-1)/RESIZE) mod(i-1,RESIZE)];
        if i == j
            graph(i,j) = 1;
            continue
        end
        dist = norm(x_i-x_j);
        if dist < r
            dist_x = exp(-dist^2/(2*sigma_x^2));
        else
            dist_x = 0;
        end
        dist_i = exp(-(img(i)-img(j))^2/(2*sigma_i^2));
        graph(i,j) = dist_i * dist_x;
    end
end

if (size(graph,1) ~= size(graph,2))
  error('Matrix not square!')
end

fprintf('End \n')
