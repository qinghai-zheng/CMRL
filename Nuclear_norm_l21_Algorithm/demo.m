clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%creating data for 2 rank  matrix
x_row=[5*ones(1,20), 4*zeros(1,2), 5*zeros(1,20), 7*zeros(1,3), 5*ones(1,5)];
x_row1=[3*zeros(1,10), 9*ones(1,2), 7*ones(1,18), ones(1,5), 9*zeros(1,15)];


X=zeros(50,50);

X(2,:) =x_row;
X(11,:) =x_row1;
X(34,:) =x_row;
X(49,:) =x_row1;
X(29,:) =x_row;
X(5,:) =x_row1;
X(21,:) =x_row;
X(44,:) =x_row1;
X(39,:) =x_row;
X(50,:) =x_row1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_image=size(X);

samp=0.4;  %%% sampling ratio 
no_meas=ceil(samp*s_image(1));  %%% calculating the  number of measurements

%%% Creating Measurement Matrix 
Mop = opBinary(no_meas,s_image(1));
Aop = Mop;
  
 for j = 2:s_image(2)
     Mop = opBinary(no_meas,s_image(1));
     Aop = opBlockDiag(1, Aop, Mop);
 end


A=opToMatrix(Aop);

y=Aop(X(:),1);

D=eye(2500);   %%% Dictionary matrix is created identity as the X is already created sparse


x_out=nuclear_l21(A,D,y,.0001,10,50,50,200); 

%%% Calculating the NMSE 
nmse = norm(X-x_out,'fro')/norm(X,'fro')
