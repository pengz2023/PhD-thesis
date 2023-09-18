%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demo for tGNBP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all') % disable all warnings
addpath('./tGNBP/')
% 1. Load data 
path = './demo_data/';
FileName_1 = [path,'GNBP_sim_p50n100_5nonzero_data(3-6)1.csv'];
xx = csvread(FileName_1);
[n,p] = size(xx);
S = xx'*xx;
FileName_2 = [path,'GNBP_sim_p50_sigmainv5(3-6).csv'];
% inverse scale matrix
CTrue = csvread(FileName_2);


% 2. Obtain the MML estimates of the hyperparameters
fprintf("MCEM algorithm for obtaining the MML estimates of the hyperparameters  - tGNBP\n");
tic;
a_lambda_sq = 0.5;
b_xi_sq = 0.5;
K = 5;
nu = 10;
EM_maxiter=10000;
[samples_a_lambda_sq, samples_b_xi_sq, samples_K, samples_nu, a_MML, b_MML, K_MML, nu_MML] = tGNBP_hyperpara_MMLcopy_withK(xx,a_lambda_sq,b_xi_sq,K,nu,EM_maxiter);
toc;
fprintf("*******************\n");


% 3. Plot the MML paths for hyperparameters
figure
l = length(samples_a_lambda_sq);
iter = 1:l;
tiledlayout(2,2);
nexttile
plot(iter,samples_a_lambda_sq);
title("a",'FontSize', 20);
nexttile
plot(iter,samples_b_xi_sq);
title("b",'FontSize', 20);
nexttile
plot(iter,samples_K);
title("K",'FontSize', 20);
nexttile
plot(iter,samples_nu); 
title("nu",'FontSize', 20);


% 4. Run the tGNBP based on MML estimates of hyperprameters
fprintf("Run tGNBP based on the MML estimates of hyperprameters\n");
tic;
a_lambda_sq = a_MML; 
b_xi_sq = b_MML; 
K = K_MML;
nu = nu_MML;
burnin=500;
nmc = 5000;
[C_save,C_vector_save,lambda_sq_save,xi_sq_save,tau_save] = tGNBP_Columnwise(xx, a_lambda_sq, b_xi_sq, K, nu, burnin, nmc);
toc;
fprintf("*******************\n");

C_save = (1-2/nu_MML)*C_save;

% 5. Use 50% credible interval for variable selection
a = 50;
C_zero = 0.5*eye(p);
for i=1:(p-1)
    for j=(i+1):p
        lb = prctile(C_save(i,j,:),(100-a)/2);
        ub = prctile(C_save(i,j,:),100-(100-a)/2);
        C_zero(i,j) = ~(lb<0 & ub>0);
    end
end
result = C_zero' + C_zero;
% Number of nonzeros in the upper/lower diagonal elements of precision matrix 
num_nonzero = (sum(result,'all')-p)/2 
% Sparsity level 
sparse_level = 1- num_nonzero/nchoosek(p,2)


% 6. Store adjacency matrix for drawing the network
path_0 = './adjacenyMatrix/';
FileName_3 = [path_0,'tGNBP_adjacenyMatrix.csv'];
writematrix(result,FileName_3);



