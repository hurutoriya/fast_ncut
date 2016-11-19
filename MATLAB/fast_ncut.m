function v = fast_ncut(A, constraint_num, MAX_ITERS, OPT_TOL)


%% Fast Normalized Cut with Linear Constraints
% Problem
% max_{v} v'Av
% Constrained Bv=c ||v|| = 1
% 
% Inputs
% - A:            Semidefinite positive matrix.
% - constraint_num: number of constrain
% - MAX_ITERS:    max number of iteration (default: 1000)
% - OPT_TOL:      optimality tolerance (default: 10^-8)
% Outputs
% - v : constraint eigen vector
% Refer to:
%   Xiao-Tong Yuan, Tong Zhang, Truncated Power Method for Sparse Eigenvalue Problems, Technical Report, 2011

fprintf('Fast Ncut ...  \n')
rng(1,'twister');

%% Init
if nargin == 0
    A = [2 0 0; 0 1 0; 0 0 0];
    MAX_ITERS = 1000;
    OPT_TOL = 10^-8;
    constraint_num = 1;
end

dim = size(A,2);

% Simple constrain: you want to same value i-th and j-th element
% You can insert 1 and -1 to i-th and j-th element.

B = zeros(constraint_num, dim-2);
B = [1 -1 B];
constraint_num = size(B,1);
c = zeros(constraint_num,1);

%% Fast ncut
k = 0;

P = eye(dim) - B'*inv(B*B')*B;
n = B'*(B*B')*c;
gamma = sqrt(1-norm(n)^2);
% if c is zero vector, v be a NaN. using any() and escape NaN.
if any(c) == 0
    v = randn(dim,1);
else
    v = gamma * P*A*n/(norm(P*A*n)) + n;
end
v = v/norm(v);
obj = v'*A*v;
obj_old = obj;
residuals = [];

tic
while k <= MAX_ITERS
    u = gamma * P*A*v/(norm(P*A*v));
    v = u + n;
    v = v/norm(v);
    obj = v'*A*v;
    % Absolute Error
    residual = abs(norm(obj - obj_old));
    residuals = [residuals residual];
    if OPT_TOL > residual
        break;
    end
    obj_old = obj;
    k = k+1;
end
toc

fprintf('Iteration: %d \n', k)
fprintf('Residual: %d \n', residual)

%% Result
% figure;
% semilogy(residuals);
% title('Residual')
% xlabel('iteration')
% ylabel('abs(v^{T}Av_{n} - v^{T}Av_{n-1})')

fprintf('Fast Ncut End \n')
end