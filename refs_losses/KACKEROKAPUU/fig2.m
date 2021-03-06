close all

Z = zeros(78,6);
j = 1;
for i=2:2:12
    Z(:,j) = Q(:,i);
    j = j+1;
end

X = zeros(78,6);
j = 1;
for i=1:2:11
    X(:,j) = Q(:,i);
    j = j+1;
end

x = zeros(78,1);
for i=1:6
    x = x + X(:,i);
end

x = x./6;



y = [40, 50, 55, 60, 65, 70];
alphas = ones(78,6);
for i = 1:78
    alphas(i,:) = y;
end


writematrix(x,'fig2_x.csv') 
writematrix(y,'fig2_y.csv') 
    
writematrix(X,'fig2_Xm.csv') 
writematrix(alphas,'fig2_Ym.csv') 
writematrix(Z,'fig2_Zm.csv') 