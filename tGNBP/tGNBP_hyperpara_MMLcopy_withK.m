% MCEM algorithm for obtaining the MML estimates of the hyperparameters in
% tGNBP

% Input:
%     Y:           n by p data matrix where n is the sample size and p is
%     the dimension;
%     a_lambda_sq: initial value of the hyperparameter 'a' in the tGNBP prior;
%     b_xi_sq:     initial value of the hyperparameter 'b' in the tGNBP prior;
%     K:           initial value of the hyperparameter 'K' in the tGNBP prior;
%     nu:          initial value of the degrees of freedom 'nu' in the tGNBP prior;
%     EM_maxiter:  maximum number of iterations for MCEM algorithm;



% Output:
%     samples_a_lambda_sq: value of hyperparameter 'a' at each iteration;
%     samples_b_xi_sq:     value of hyperparameter 'b' at each iteration;
%     samples_K:           value of hyperparameter 'K' at each iteration;
%     samples_nu:          value of degrees of freedom 'nu' at each iteration;
%     a_MML:               the maximum marginal likelihood (MML) estimate of 'a';
%     b_MML:               the maximum marginal likelihood (MML) estimate of 'b';
%     K_MML:               the maximum marginal likelihood (MML) estimate of 'K';
%     nu_MML:              the maximum marginal likelihood (MML) estimate of 'nu';


function [samples_a_lambda_sq, samples_b_xi_sq, samples_K, samples_nu, a_MML, b_MML, K_MML, nu_MML] = tGNBP_hyperpara_MMLcopy_withK(Y,a_lambda_sq,b_xi_sq,K,nu,EM_maxiter) 
[n,p] = size(Y);
lambda_sq_save = zeros(p,p,EM_maxiter);
xi_sq_save = zeros(p,p,EM_maxiter);
tau_save = zeros(n,EM_maxiter);

% Initialization 
Sig = eye(p); C = eye(p);
lambda_sq(1:p,1:p)=1; xi_sq(1:p,1:p)=1;
tau(1:n,1) = 1;

% Maximum number of iterations for EM algorithm
%     EM_maxiter = burnin;
%     EM_tol = 1e-04; % Tolerance for the EM algorithm
    EM_tol = 1e-03; % Tolerance for the EM algorithm
    EM_dif = 1;     % Initialize
    abKnu_current = [a_lambda_sq, b_xi_sq, K,nu];
    abKnu_update = [a_lambda_sq, b_xi_sq, K,nu];
    samples_a_lambda_sq = [a_lambda_sq];
    samples_b_xi_sq = [b_xi_sq];
    samples_nu = [nu];
    samples_K = [K];
    
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



for iter = 1:EM_maxiter
            
       if(mod(iter,1000)==0)
        fprintf('iter = %d \n',iter);
       end
    
    % S_tau
    S = Y' * diag(tau) * Y;
    
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
        
        %% Below updating Covariance matrix according to one-column change of precision matrix
        invC11beta = invC11*beta;
        Sig(ind_noi,ind_noi) = invC11+invC11beta*invC11beta'/gam;
        Sig12 = -invC11beta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;
        
    % Step 2. Sample lambda_sq off-diagonal
    Cadjust = max(abs(beta),10^-8);
    xi_sq_adjust = max(abs(xi_sq(ind_noi,i)),10^-8);
    ppost = -1/2 + a_lambda_sq;
    apost = 2;
    bpost =  (Cadjust.^2)./xi_sq_adjust;   
    lambda_sq(ind_noi,i) = gigrndHandle(ppost,apost,bpost);
    lambda_sq(i,ind_noi) = lambda_sq(ind_noi,i)';
    
    % Step 3. sample xi_sq off-diagonal 
    shape = 1/2 + b_xi_sq;
    scale = lambda_sq(ind_noi,i)./(0.5*Cadjust.^2+lambda_sq(ind_noi,i));
    xi_sq(ind_noi,i) = 1./gamrnd(shape, scale);
    xi_sq(i,ind_noi) = xi_sq(ind_noi,i)';
    
    end
    
%%% Step 4. sample tau
    % rate=(diag(Y*C*Y')+nu)./2
    scale_tau = 2./(diag(Y*C*Y')+nu); 
    tau = gamrnd((p+nu)/2, scale_tau);
    
    % store values 
     C_save(:,:,iter) = C;
     lambda_sq_save(:,:,iter) = lambda_sq;
     xi_sq_save(:,:,iter) = xi_sq;
     tau_save(:,iter) = tau;

    % EM Monte Carlo to find optimal a_lambda_sq, b_xi_sq, K and nu
    if mod(iter,100)==0
        abKnu_current = [a_lambda_sq, b_xi_sq, K,nu]; % Previous ab update gets saved here
        low = iter-99;
        high = iter;
        
%       update  a_lambda_sq
        ln_lambda_terms = mean(log(lambda_sq_save(:,:,low:high)),3);
        fa = @(a) -p*(p-1)*psi(a) + sum(ln_lambda_terms,"All");
        a_lambda_sq = fzero(fa,[eps,1e16]);
        samples_a_lambda_sq(iter/100+1) = a_lambda_sq;
  
%         % update  b_xi_sq
        ln_xi_terms = mean(log(xi_sq_save(:,:,low:high)),3);
        fb = @(b) p*(p-1)*psi(b) + sum(ln_xi_terms,"All");
        b_xi_sq = fzero(fb,[eps,1e16]);
        samples_b_xi_sq(iter/100+1) = b_xi_sq;  
        
        
        % update K
        K = 2.0*p/sum(diag(mean(C_save(:,:,low:high),3)));
        samples_K(iter/100+1) = K; 
        
        % update nu
        ln_tau_terms = mean(log(tau_save(:,low:high)),2);
        tau_terms = mean(tau_save(:,low:high),2);
        fnu_half = @(nu_half) log(nu_half) - psi(nu_half) - mean(tau_terms-ln_tau_terms,'all') + 1;
        nu = max(2*fzero(fnu_half,[eps,1e16]),2.1); %ensure nu>2
        samples_nu(iter/100+1) = nu;  
        
        
        % calculate difference
        abKnu_update = [a_lambda_sq, b_xi_sq, K,nu]; 
        diffe = abKnu_update-abKnu_current;
        
        % update the difference
        EM_dif = sqrt(sum(diffe.^2));
        
        if EM_dif<EM_tol
            disp(['******Excute Break****** ',num2str(EM_dif),' at the iteration: ',num2str(iter), ' >>>>>>']);
            break;
        end
    end
end
            a_MML = a_lambda_sq;
            b_MML = b_xi_sq;
            K_MML = K;
            nu_MML = nu;
            disp(['******a_MML: ',num2str(a_MML),' ***']);
            disp(['******b_MML: ',num2str(b_MML),' ***']);
            disp(['******K_MML: ',num2str(K_MML),' ***']);
            disp(['******nu_MML: ',num2str(nu_MML),' ***']);
end

