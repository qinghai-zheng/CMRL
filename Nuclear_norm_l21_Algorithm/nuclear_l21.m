function x_out= nuclear_l21(A,D,Y,lambda1,lambda2,M,N,iter_count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function description:
%output argument x_out: recovered signal
%input arguments    A: measurement matrix
%                   D: Dictionary matrix which is taken to identity in case
%                      of synthetic data otherwise transform domain in which
%                      the signal is sparse
%                   Y: Observed signal
%             lambda1: tuning paramter for nuclear norm term 
%             lambda2: tuning parameter for l21 term
%             M,N     :Dimensions of the signal MXN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Statement:
%    lambda1||x||* +lambda2||Dx||2,1 subjected to y=Ax

% This function solves for the recovery of a signal which is low rank and
% joint sparse using split Bregman mathod.The above problem takes the form
% of 1/2 ||y-Ax||^2 +lambda1||W||* +lambda2||DZ||2,1 +eta1/2||W-X-B1||^2
% +eta2/2||Z-X-B2||^2
% B1 and B2 are Bregman variable
% W and Z are the proxy variables involved in splitting
% eta1 and eta2 are regularization parameter

 
%%% values of interanl tuning parameter
eta1=.0001;
eta2=.0001;

%%%% Desired tolerence for signal recovery (One can change this value accordingly)
tol=1e-5; 


%%% Caculating the size from observed signal  for sizes of
%%% proxy  variables Wand Z and bregman variable B1 and B2
s_y=size(Y);
y_vec=Y;
sa2=M*N;

%%%%Initialization of the proxy  and bregman variable invoved in Split Bregman 
W=zeros(sa2,1);
Z=zeros(sa2,1);
B1=ones(sa2,1);
B2=ones(sa2,1);


%%%%%% Values of Algorithm Internal Parameters
alpha=1.01*max(eig(A'*A));
c=1.05*max(eig(D'*D));
%%%%%%%%%%%%%%%%%%%


%%%%%%%Creating a variable for the values of objective function at different iterations
objective_func=[]; 


%%%%% Algorithm steps %%%%%%%%%
 for iteration=1:iter_count
     %%% step for solving ||y-Ax||^2%%%%%%%
     X=l2_term(A,W,Z,y_vec,B1,B2,eta1,eta2); 
     
     %%%% Updating W%%%%%%%
     W=X+B1;
     
     %%% Step for nuclear norm term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [W,s_thld]=nuclear_norm(W,sa2,M,N,lambda1,eta1);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %%% Step for l21 term taking care of joint sparsity in the data%%%%%
     [Z,l2_row] = l21_term(X,Z,D,B2,alpha,lambda2,sa2,M,N,c);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %%% Update for Bregman Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     B1 = B1+X-W;
     B2=B2+X-Z;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     objective_func(iteration)=norm(y_vec-A*X)+lambda1*s_thld+lambda2*sum(abs(1./l2_row));
     
     if norm(y_vec-A*X)< tol  %%% exit criteria if desired tolerance is achieved 
          break;
     end
      
 end
 
  x_out = reshape(X,M,N); %%% recovered signal 
  
   %%% plot for objective function
   plot(objective_func)  
   title('Plot of Objective Function')
   xlabel('No of iterations')
   ylabel('Value of Objective Function')
end 












   