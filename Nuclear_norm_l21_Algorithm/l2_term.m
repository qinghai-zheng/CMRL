function x_temp=l2_term(A,W,Z,Y,B1,B2,eta1,eta2)
%%%%% x_temp=l2_term(A,W,Z,Y,B1,B2,eta1,eta2) Function to solve least
%%%%% square problem analytically

temp=(A'*A+eta1*eye(size(A,2)) + eta2*eye(size(A,2)));
temp1=eta1*(W-B1)+eta2*(Z-B2)+A'*Y;
x_temp=temp\temp1;

end