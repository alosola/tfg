close all

Q = table2array(readtable('staggerangle.csv'));

Z = zeros(43,7);
j = 1;
for i=2:2:14
    Z(:,j) = Q(3:end,i);
    j = j+1;
end
Z = flip(Z,2)
Z = flip(Z,1)



X = zeros(43,7);
j = 1;
for i=1:2:14
    X(:,j) = Q(3:end,i);
    j = j+1;
end
X = flip(X,2)

x = zeros(43,1);
for i=1:6
    x = x + X(:,i);
end

x = x./6;



y = 50:5:80;
alphas = ones(43,7);
for i = 1:43
    alphas(i,:) = y;
end


writematrix(x,'stg_x.csv') 
writematrix(y,'stg_y.csv') 
    
writematrix(X,'stg_Xm.csv') 
writematrix(alphas,'stg_Ym.csv') 
writematrix(Z,'stg_Zm.csv') 