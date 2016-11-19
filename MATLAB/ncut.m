function clusters = ncut(L, D, k)
%% Normalized Cut
% Inputs
% - L : Laplacian Matrix.
% - D : Diagonal Matrix.
% - k : number of cluster.
% Output
% - clusters : spectral clustering result.
% Refer to:
% Normalized Cuts and Image Segmentation Jianbo Shi and Jitendra Malik IEEE Transactions on Pattern Analysis and Machine Intelligence(PAMI) 2000
% Sample Code link : http://www.cis.upenn.edu/~jshi/software/

%% Normalized Cuts
fprintf('Normalized Cut ...  \n')

% Computing Smallest 5 eigen vectors

[V,D] = eigs(L,D,5,'sa');

% se : second smallest eigenvectors.
se = V(:,2:3);
clusters = kmeans(se,k);

fprintf('End \n')

end