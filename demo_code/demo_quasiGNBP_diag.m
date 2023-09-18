%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Demo for quasiGNBP-diag (quasiGNBP algorithm when the diagonal elements are known) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all') % disable all warnings
addpath('./quasiGNBP_diag/')
% 1. Load data  
path = './demo_data/';
FileName_1 = [path,'GNBP_sim_p50n100_5nonzero_data(3-6)1.csv'];
xx = csvread(FileName_1);
[n,p] = size(xx);
S = xx'*xx;
FileName_2 = [path,'GNBP_sim_p50_sigmainv5(3-6).csv'];
CTrue = csvread(FileName_2);


% 2. Obtain the MMQL estimates of the hyperparameters
fprintf("MMQL Search Hyperparameters quasiGNBP with unknown diagonal elements\n");
tic;
a_lambda_sq = 0.5;
b_xi_sq = 0.5;
EM_maxiter=30000;
% Use the diagonal elements of true precision matrix
Diag = diag(CTrue);  
[samples_a_lambda_sq, samples_b_xi_sq, a_MML, b_MML] = quasiGNBP_Diag_hyperpara_MML_withoutK(S,Diag,a_lambda_sq,b_xi_sq,EM_maxiter);
toc;
fprintf("*******************\n");


% 3. Plot the MMQL paths for hyperparameters
figure
l = length(samples_a_lambda_sq);
iter = 1:l;
tiledlayout(1,2);
nexttile
plot(iter,samples_a_lambda_sq);
title("a",'FontSize', 20);
nexttile
plot(iter,samples_b_xi_sq);
title("b",'FontSize', 20);


% 4. Run the quasiGNBP-diag based on MMQL estimates of hyperprameters
fprintf("Run quasiGNBP-diag algorithm based on the MMQL estimates of hyperprameters\n");
tic;
a_lambda_sq = a_MML; 
b_xi_sq = b_MML; 
burnin=500;
nmc = 5000;
[quasiGNBP_C_save,lambda_sq_save,xi_sq_save] = quasiGNBP_Diag_elementwise_withoutK(S, Diag, a_MML, b_MML, burnin, nmc);
toc;
fprintf("*******************\n");

% 5. Use 50% credible interval for variable selection
a = 50;
percen_lb = (100-a)/2;
percen_ub = 100-(100-a)/2;
quasiGNBP_C_percen = prctile(quasiGNBP_C_save, [percen_lb percen_ub], 3);
% Note that zero pattern matrix isZero could be asymmetric; 
isZero = quasiGNBP_C_percen(:,:,1)<0 & quasiGNBP_C_percen(:,:,2)>0; 
% If the 50% posterior credible interval for C(i,j) or C(j,i) contains zero, the corresponding 
% element is considered to be zero. In this way, the reuslting zero pattern of C could be asymmetric; 
% According to Zhang et al. (2022), C(i,j) is considered as nonzero if either C(i,j) or C(j,i) 
% is nonzero. In other words, C(i,j) is considered as zero iff both C(i,j) and C(j,i) is zero. 
quasiGNBP_C_NonZero = ~((isZero + isZero')==2);
diag(quasiGNBP_C_NonZero) = 1;
result = quasiGNBP_C_NonZero;
% Number of nonzeros in the upper/lower diagonal elements of precision matrix 
num_nonzero = (sum(result,'all')-p)/2 
% Sparsity level
sparse_level = 1- num_nonzero/nchoosek(p,2)



% 6. Store adjacency matrix for drawing the network
path_0 = './adjacenyMatrix/';
FileName_3 = [path_0,'quasiGNBP_diag_adjacenyMatrix.csv'];
writematrix(result,FileName_3);

% References:
% Ruoyang Zhang, Yisha Yao, and Malay Ghosh. Contraction of a quasi-bayesian model with shrinkage priors in precision matrix estimation. Journal of Statistical Planning and Inference, 221:154–171, 2022.