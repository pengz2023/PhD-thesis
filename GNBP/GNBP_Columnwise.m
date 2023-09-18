% Gibbs sampler for GNBP


% Input:
%     S = Y'*Y :   sample covariance matrix * n;
%     n:           sample size;
%     a_lambda_sq: hyperparameter 'a' in the GNBP prior;
%     b_xi_sq:     hyperparameter 'b' in the GNBP prior;
%     K:           hyperparameter 'K' in the GNBP prior;
%     burnin:      number of MCMC burnins;
%     nmc :        number of saved MCMC samples;
%     

% Output:
%     C_save: p by p by nmc matrices of saved posterior samples of
%     precision matrix;

%     C_vector_save: p*(p-1)/2 by nmc vectors of saved posterior samples of
%     upper diagonal elements of precision matrix;

%     lambda_sq_save: p by p by nmc matrices of saved samples of lambda
%     squared (latent variables);

%     xi_sq_save: p by p by nmc matrices of saved samples of xi
%     squared (latent variables);


function [C_save,C_vector_save,lambda_sq_save,xi_sq_save] = GNBP_Columnwise(S,n,a_lambda_sq,b_xi_sq,K,burnin,nmc) 

[p] = size(S,1); C_save = zeros(p,p,nmc); 
C_vector = zeros(p*(p-1)/2,1);
C_vector_save = zeros(p*(p-1)/2,nmc);
lambda_sq_save = zeros(p,p,nmc);
xi_sq_save = zeros(p,p,nmc);


% Initialization 
Sig = eye(p); C = eye(p);
lambda_sq(1:p,1:p)=1; xi_sq(1:p,1:p)=1;

ind_noi_all = zeros(p-1,p);
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       
       ind_noi_all(:,i) = ind_noi;
end



for iter = 1: burnin+nmc    
            
       if(mod(iter,1000)==0)
        fprintf('iter = %d \n',iter);
       end

    %%% Step 1. sample Sig and C = inv(Sig)        
    for i = 1:p
        ind_noi = ind_noi_all(:,i);
        lambda_sq_temp = lambda_sq(ind_noi,i);
        xi_sq_temp = xi_sq(ind_noi,i);
        tau_temp = lambda_sq_temp.* xi_sq_temp;
        Sig11 = Sig(ind_noi,ind_noi); Sig12 = Sig(ind_noi,i);
        invC11 = Sig11 - Sig12*Sig12'/Sig(i,i);
        Ci = (S(i,i)+K)*invC11+diag(1./tau_temp);
        Ci_chol = chol(Ci);
        mu_i = -Ci\S(ind_noi,i);
        beta = mu_i+ Ci_chol\randn(p-1,1);
        C(ind_noi,i) = beta;
        C(i,ind_noi) = beta;
        gam = gamrnd(n/2+1,2/(S(i,i)+K)); % rate=(s_22+k)/2
        C(i,i) = gam+beta'*invC11*beta;
        %%%%
  
        %% Below updating Covariance matrix according to one-column change of precision matrix
        invC11beta = invC11*beta;
        Sig(ind_noi,ind_noi) = invC11+invC11beta*invC11beta'/gam;
        Sig12 = -invC11beta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;
 % %%% Step 2. Sample lambda_sq off-diagonal
        Cadjust = max(abs(beta),10^-7);
        xi_sq_adjust = max(abs(xi_sq(ind_noi,i)),10^-7);
        ppost = -1/2 + a_lambda_sq;
        apost = 2;
        bpost =  (Cadjust.^2)./xi_sq_adjust;   
        lambda_sq(ind_noi,i) = gigrndHandle(ppost,apost,bpost);
        lambda_sq(i,ind_noi) = lambda_sq(ind_noi,i)';

%%% Step 3. sample xi_sq off-diagonal 
        shape = 1/2 + b_xi_sq;
        scale = lambda_sq(ind_noi,i)./(0.5*Cadjust.^2+lambda_sq(ind_noi,i));
        xi_sq(ind_noi,i) = 1./gamrnd(shape, scale);
        xi_sq(i,ind_noi) = xi_sq(ind_noi,i)';
    
   end
       if iter >burnin
            C_save(:,:,iter-burnin) = C;
            C_vector = C(tril(true(size(C)),-1));
            C_vector_save(:,iter-burnin) = C_vector;
            lambda_sq_save(:,:,iter-burnin) = lambda_sq;
            xi_sq_save(:,:,iter-burnin) = xi_sq;
       end
end
end

