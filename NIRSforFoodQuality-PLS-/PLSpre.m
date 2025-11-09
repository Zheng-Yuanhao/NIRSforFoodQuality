load('PLS20221109.mat')
load('Redcitrus_650-1050_45_20221110_val.mat')
X=X(:,1:1800);
m=size(X,1);
regcoef_original=PLS.regcoef_original;
n=size(regcoef_original,1);
Y_val=zeros(m,1);
for i=1:m
    y_val=0;
    for j=1:n
        if j<n
            y_val=y_val+X(i,j)*regcoef_original(j);
        else
            y_val=y_val+regcoef_original(j);
        end
    end 
    Y_val(i)=y_val;
end