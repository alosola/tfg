
P = table2array(readtable('dataset_machhub.csv'));
M2H = P(:,1);

fine = 500;

Av = linspace(0.01,20,fine);
Bv = linspace(0.01,20,fine);
Cv = linspace(0.01,20,fine);


for i=1:fine
    A = Av(i)
    
    for j = 1:fine
        B = Bv(i);
        
        for k = 1:fine
            C = Cv(i);
            
            dp = find_dp(A,B,C,M2H);
            
            diff = dp - P(:,2);
            
            if max(diff) < 1e-3
                solution = diff;
                break
            end
        end
    end
end
                

figure(1)
plot(P(:,1),P(:,2),'.k')
xlim([0.2,1.2])
ylim([0,5])

function dp = find_dp(A, B, C, M2H)
dp = A.*(M2H-B).^C;

end