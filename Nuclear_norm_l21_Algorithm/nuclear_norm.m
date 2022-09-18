function [w_temp,sumS]= nuclear_norm(W,sa2,M,N,lambda1,eta1)

%%%%% function for minimizing nuclear norm%%%%
     W_mat=reshape(W,M,N);
     [U,S,V]=svd(W_mat);
     sumS=0;
      for i=1:min(size(S))
         S(i,i) = soft(S(i,i),(4*lambda1)/eta1);
         sumS=sumS+S(i,i);
      end
      W_mat=U*S*V';

     w_temp=reshape(W_mat,sa2,1);
end