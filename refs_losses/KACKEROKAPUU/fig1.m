close all


P = zeros(75, 2, 7);
P(:,:,1) = table2array(readtable('fig1_alpha80.csv'));
P(:,:,2) = table2array(readtable('fig1_alpha75.csv'));
P(:,:,3) = table2array(readtable('fig1_alpha70.csv'));
P(:,:,4) = table2array(readtable('fig1_alpha65.csv'));
P(:,:,5) = table2array(readtable('fig1_alpha60.csv'));
P(:,:,6) = table2array(readtable('fig1_alpha50.csv'));
P(:,:,7) = table2array(readtable('fig1_alpha40.csv'));

for i=1:7
    figure(1)
    scatter(P(:,1,i),P(:,2,i),'.');
    hold on
    
end


%%


alphas = ones(75,3);
for i = 1:3
    alphas(:,i) = (30+i*10);
end
for i = 4:7
    alphas(:,i) = (45+i*5);
    
end

B = zeros(75, 3, 7);
for i=1:7
    B(:,1,i) = alphas(:,i);
    B(:,2:3,i) = P(:,:,i);
    
end

B(75,2,7) = B(74,2,7)

% for i = 1:7
%     B(:,2,i) = B(:,2,1)
% end

A = zeros(75*7,3)
aux = B(:,:,1)
for i=2:7
    A = [aux;B(:,:,i)];
    aux = A
end

for i=1:7
    figure(1)
    scatter(B(:,2,i),B(:,3,i),'.');
    hold on
    
end

figure(2)
scatter3(A(:,1,:),A(:,2,:),A(:,3,:))

figure(3)
hold on
aux = 0
for i=1:7
    plot(B(:,2,i),'.')
    aux = aux + B(:,2,i)
end

x = (aux/7)'
plot(x,'*')

%%
% x = B(:,2,i)';
y = flip(alphas(1,:));

z = zeros(7,75)
for i =1:7
    z(i,:) = B(:,3,i);
end

answer =  interp2(x, y, z, 0.6, 35, 'spline');
   
    
figure(4)
surf(x,y,z)

writematrix(x,'fig1_x.csv') 
writematrix(y,'fig1_y.csv') 
    

%%
X = zeros(75,7)
for i=1:i
    X(:,i) = B(:,2,i)
end

writematrix(X,'fig1_Xm.csv') 
writematrix(alphas,'fig1_Ym.csv') 
writematrix(z','fig1_Zm.csv') 
    
    
    
    
    