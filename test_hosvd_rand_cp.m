clc;
clear;
close all;

addpath('.../tensor_toolbox-v3.6')%replace ... with path to the toolbox

load('tmax.mat')
load('tmin.mat')
load('pr_trans.mat')

X = cat(3, pr_trans, tmin, tmax);  % Stack matrices along 3rd dimension
T = tensor(X);

% Run CPD-ALS
[M, output] = cp_als(T, 3,'init','nvecs', 'printitn',1); %for HOSVD
%[M, output] = cp_als(T, 3,'init','random', 'printitn',1); %for random

% Save results
%save('cp_als_3.mat', 'M', 'output');

% Build estimated model from normalized updated matrices
A = M.u{1};   % size: size(X,1) x R
B = M.u{2};   % size: size(X,2) x R
C = M.u{3}; 
Co = M.lambda;

    X_est = tensor_create2(A,B,C,Co);
    relative_cost = norm_fro(X-X_est)/norm_fro(X)

