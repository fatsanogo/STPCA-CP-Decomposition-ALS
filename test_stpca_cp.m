%Fatou Sanogo created this file in Aug 2025 

clc;
clear;
close all;

addpath('.../tensorlab_2016-03-28') %replace ... with path to the toolbox

load('tmax.mat')
load('tmin.mat')
load('pr_trans.mat')
load ('stpca_init3.mat'); %for rank 2, replace 2 with 3 for rank 3

X = cat(3, pr_trans, tmin, tmax);  % Stack matrices along 3rd dimension

% Initialize factor matrices
U0 = {U_time, U_space, U_var};

% Format the tensor 
T = fmt(X, true); 

% Run CPD-ALS
[U, output] = cpd_als(T, U0);

res2_final   = 2*output.fval(end);        % ||T - T_hat||_F^2 at the last iter
T2           = frob(T,'squared');         % ||T||_F^2  (consistent for structured T)
relerr_final = sqrt(res2_final / T2);     % final relative error (fraction)
fit_final_pct = 100*(1 - relerr_final);   % final fit in percent

fprintf('Final rel. error = %.6f (%.2f%%), fit = %.2f%%\n', ...
        relerr_final, 100*relerr_final, fit_final_pct);

% Save results
%save('cpd_als_30.mat', 'U', 'output');