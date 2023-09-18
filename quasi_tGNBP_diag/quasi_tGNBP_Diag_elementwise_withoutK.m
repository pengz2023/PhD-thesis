% Input:
%     Y:           n by p data matrix where n is the sample size and p is
%     the dimension;
%     Diag:        diagonal elements of the inverse scale matrix;
%     a_lambda_sq: hyperparameter 'a' in the quasi-tGNBP prior;
%     b_xi_sq:     hyperparameter 'b' in the quasi-tGNBP prior;
%     nu:          degrees of freedom;
%     burnin:      number of MCMC burnins;
%     nmc :        number of saved MCMC samples;


% Output:
%     C_save: p by p by nmc matrices of saved posterior samples of
%     inverse scale matrix;

%     lambda_sq_save: p by p by nmc matrices of saved samples of lambda
%     squared (latent variables);

%     xi_sq_save: p by p by nmc matrices of saved samples of xi
%     squared (latent variables);

%     tau_save:   n by nmc vectors of saved samples of tau (latent variables);



function [C_save,lambda_sq_save,xi_sq_save,tau_save] = quasi_tGNBP_Diag_elementwise_withoutK(Y,Diag, a_lambda_sq,b_xi_sq,nu,burnin,nmc) 

[n,p] = size(Y);
C_save = zeros(p,p,nmc); 
lambda_sq_save = zeros(p,p,nmc);
xi_sq_save = zeros(p,p,nmc);
tau_save = zeros(n,nmc); 

% Initialization 
C = diag(Diag);
lambda_sq(1:p,1:p)=1; xi_sq(1:p,1:p)=1;
index = true(p);
index(1:(p+1):end)= 0;
tau(1:n,1) = 1;

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


for iter = 1:burnin+nmc
    if(mod(iter,1000)==0)
        fprintf('iter = %d \n',iter);
    end
% 
    
    % setting S_tau
    S = Y' * diag(tau) * Y;
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


    %%%  sample tau
    central_mat = C * diag(1./(diag(C))) * C';
     % rate=(diag(Y*central_mat*Y')+nu)./2
    scale_tau = 2./(diag(Y*central_mat*Y') + nu); 
    tau = gamrnd((p+nu)/2, scale_tau);
    
    if iter >burnin
       C_save(:,:,iter-burnin) = C;
       lambda_sq_save(:,:,iter-burnin) = lambda_sq;
       xi_sq_save(:,:,iter-burnin) = xi_sq;
       tau_save(:,iter-burnin) = tau;
    end

end

end

