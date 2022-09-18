function vaule = Norm21(A)

[x,y] = size(A);
value_tmp = 0;
vaule = 0;
for i = 1:y
    for j = 1:x
        value_tmp = value_tmp + A(j,i)*A(j,i);
    end
    vaule = vaule + value_tmp^(0.5);
    value_tmp = 0; 
end


