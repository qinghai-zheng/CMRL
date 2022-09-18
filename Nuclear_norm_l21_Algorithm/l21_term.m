function [z_out,sig]=l21_term(X,Z,D,B2,alpha,lambda2,sa2,M,N,c)
%%%%function [z_out,sig]=l21_term(X,Z,D,B2,alpha,lambda2,sa2,M,N,c)
% function for solving L21 part of the problem
sig=[];
z_temp=X+B2;

r=zeros(size(Z));
z_in=zeros(size(B2,1),1);
 
       z_in=z_temp;
       b=z_in;
       z_mat=reshape(z_temp,M,N);
       T=z_mat;
                for in=1:size(T,1)
                    sig(in,:)=(1/norm(T(in,:)));
                end

      temp=diag(1./((alpha/lambda2)*(1./sig)+c));
      t1=(c*r+D*(b-D'*r));
      ret1=reshape(t1,M,N);
      r=temp*ret1;
      r=reshape(r,sa2,1);
      z_temp=b-D'*r;

  z_out=z_temp;
end