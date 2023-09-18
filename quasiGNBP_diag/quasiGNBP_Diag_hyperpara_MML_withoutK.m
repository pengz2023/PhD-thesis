% Input:
%     S = Y'*Y :   sample covariance matrix * n;
%     Diag:        diagonal elements of the precision matrix;
%     a_lambda_sq: initial value of the hyperparameter 'a' in the quasiGNBP prior;
%     b_xi_sq:     initial value of the hyperparameter 'b' in the quasiGNBP prior;
%     EM_maxiter:  maximum number of iterations for MCEM algorithm;


% Output:
%     samples_a_lambda_sq: value of hyperparameter 'a' at each iteration;
%     samples_b_xi_sq:     value of hyperparameter 'b' at each iteration;
%     a_MML:               the maximum marginal quasi-likelihood (MMQL) estimate of 'a';
%     b_MML:               the maximum marginal quasi-likelihood (MMQL) estimate of 'b';


function [samples_a_lambda_sq, samples_b_xi_sq, a_MML, b_MML] = quasiGNBP_Diag_hyperpara_MML_withoutK(S,Diag,a_lambda_sq, b_xi_sq, EM_maxiter) 

[p] = size(S,1);  
lambda_sq_save = zeros(p,p,EM_maxiter);
xi_sq_save = zeros(p,p,EM_maxiter);


% Initialization 
% Sig = eye(p); 
C = diag(Diag);
lambda_sq(1:p,1:p)=1; xi_sq(1:p,1:p)=1;
index = true(p);
index(1:(p+1):end)= 0;


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

% Maximum number of iterations for EM algorithm
%     EM_maxiter = burnin;
%     EM_tol = 1e-06; % Tolerance for the EM algorithm
    EM_tol = 1e-03; % Tolerance for the EM algorithm
    EM_dif = 1;     % Initialize
    ab_current = [a_lambda_sq, b_xi_sq];
    ab_update = [a_lambda_sq, b_xi_sq];
    samples_a_lambda_sq = [a_lambda_sq];
    samples_b_xi_sq = [b_xi_sq]; 


for iter = 1:EM_maxiter
       if(mod(iter,1000)==0)
        fprintf('iter = %d \n',iter);
       end

    for i = 1:p
        for j = ind_noi_all(:,i)'            
            num = C(:,i)'*S(:,j)- C(j,i)*S(j,j);
            denom = S(j,j) + C(i,i)/(lambda_sq(j,i)*xi_sq(j,i));
            C(j,i) = -num/denom + sqrt(C(i,i)/denom)*randn(1);
        end
    end

    %   Sample lambda_sq off-diagonal
    Cadjust = max(abs(C(index)),10^-7);
    xi_sq_adjust = max(abs(xi_sq(index)),10^-7);
    ppost = -1/2 + a_lambda_sq;
    bpost =  (Cadjust.^2)./xi_sq_adjust;
    lambda_sq(index) = gigrndHandle(ppost,2,bpost);

    %%% Sample xi_sq off-diagonal 
    shape = 1/2 + b_xi_sq;
    scale = lambda_sq(index)./(0.5*Cadjust.^2+lambda_sq(index));
    xi_sq(index) = 1./gamrnd(shape, scale);
    
    % store values
     lambda_sq_save(:,:,iter) = lambda_sq;
     xi_sq_save(:,:,iter) = xi_sq;

    % EM Monte Carlo to find optimal a_lambda_sq, b_xi_sq and K
    if mod(iter,100)==0
        ab_current = [a_lambda_sq, b_xi_sq]; % Previous ab update gets saved here
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
        
        
        % calculate difference
        ab_update = [a_lambda_sq, b_xi_sq]; 
        diffe = ab_update-ab_current;
                
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
            disp(['******a_MML: ',num2str(a_MML),' ***']);
            disp(['******b_MML: ',num2str(b_MML),' ***']);
end

