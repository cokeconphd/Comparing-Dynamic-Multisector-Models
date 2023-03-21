%Figure 1: MPY17: A Multisector DSGE Model with Intratemporal Adjustment
%Costs of Capital

rho1=-1;
rho2=-2;
rho3=-3;
K=1;
k1=0.00:0.01:1;
k21=K-k1.^(-rho1);
k22=K-k1.^(-rho2);
k23=K-k1.^(-rho3);

h=figure;
plot(k1,k21,'LineWidth' , 2)
hold on
plot(k1,k22,'LineWidth' , 2)
hold on
plot(k1,k23,'LineWidth' , 2)
hold off
%title('Isoquant K=1')
xlabel('k_1')
ylabel('k_2')
legend('\rho=-1','\rho=-2','\rho=-3')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\intratemp_adj.png');

