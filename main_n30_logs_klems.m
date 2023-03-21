%Miranda-Pinto and Young 2019: Solves a N sector 1 capital model DSGE model
%with sectoral linkages, calibrated for the U.S. sectors (KLEMS)

clearvars
clc  
tic

SOLVE=1;%1 is need to solve intratemporal adj cost model
CORRELATIONS=0; %1 to calculate implied distribution of pairwise correlation across sectors
VAR=0;% 1 to estimate VAR(1) data and models and then IRFs
GRAPHS=0    ; %1 to plot the VAR(1) graphs
SPECTRALST=0; %1 to estimate barltett density for stationary TFP

% Parameter Values

BETTA = 0.96;     %discount rate
DELTA = 0.1;    %depreciation rate (annual)
SIGM  = 1;       %CRRA coeficient Utility function
RHO   = -1.1              ;       %Intratemporal adjustment cost of capital parameter (-1 means no cost)
VARRHO= 0.99;     %persistence of productivity shock
SIGZ  =0.0076;    %volatility of techmnology shock (this will have to be a matriz later on
PSSI  = 1;        %Coeficient for labor/leiusure utility

%Load the IO shares and production function shares

 %av_theta=csvread('avtheta.csv');
 %av_theta=av_theta+10e-8;
 %save ('av_theta.mat','av_theta');

load av_theta %av_theta KLEMS Atalay 07 AEJ: Macro
av_theta_row=av_theta;
load gamma_klems %I-O matrix Atalay 07 AEJ: Macro
GAMMA=gamma_M_data;
load alpha_klems %capital shares Klems (labor is 1-alpha)
ALFA=alpha;
load mu_klems.mat
GAMMA=(mu_vector_data*(ones(1,30)))'.*GAMMA;
ioshares=sum(GAMMA)';

[N, col]=size(GAMMA);

for i=1:col
    for j=1:col
        if GAMMA(i,j)==0
    GAMMA(i,j) = 10e-20;
        end
    end
end

%Other parameters for simulations
shry=csvread('output_shares_klems.csv');
 
T=500; %size of the simulated sample from models
SIM=1; %number of simulation for each VAR estimation
w=0:0.01:2.14; %frequencies
nimp=1000; %number of period to estimate population covariance matrix Y

% Steady State and Parameter/Matrix Values

VARFI=ones(1,N);
Mat_mlamda = kron(ones(N,1), eye(N)) - kron(eye(N),ones(N,1)); 
my=kron(ones(N,1),eye(N));
alphad = diag(alpha,0);
phi = eye(N) - alphad - repmat(sum(GAMMA)', 1, N).*eye(N);
gma_tilde0 = (kron(GAMMA',ones(1,N)));
gma_tilde1 = kron(ones(1,N),eye(N));
gma_tilde  =  gma_tilde0.*gma_tilde1;
z=ones(N,1);
alrho=alphad*(RHO+1)/RHO;
C_l = phi*ones(N,1);
C_beta = ((1-BETTA*(1-DELTA))./(BETTA*ALFA)); 

I=eye(N);

BETTA_TILDE=1-BETTA*(1-DELTA);
M_c=SIGM*(kron(eye(N),ones(N,1)) - kron(ones(N,1),eye(N))); 
S_m=kron(eye(N),ones(1,N)); 

%Steady State Special Case rho=-1

temp = alphad*(I-ones(N,1)*av_theta_row')+gma_tilde*Mat_mlamda + phi;  
temp = inv(temp);
ln_lamda = -temp*(z+alphad*(ones(N,1)*av_theta_row'*log(av_theta_row)-log(C_beta))+gma_tilde*reshape(log(GAMMA'),1,N*N)' - phi*(log(PSSI)*ones(N,1)) + phi*log(C_l) );
LAM = exp(ln_lamda);
ln_mu=av_theta_row'*ln_lamda-av_theta_row'*log(av_theta_row);
M_XI=exp(ln_mu)*av_theta_row./LAM;
LAM301=LAM;
save ('LAM301.mat','LAM301');
CC=LAM.^(1/-SIGM);

lamda_rep=repmat(LAM,1,N)'./repmat(LAM,1,N);
Mat_my = (GAMMA.*lamda_rep); 

C_exp=diag(exp((I-ones(N,1)*av_theta_row')*log(LAM)+ones(N,1)*av_theta_row'*log(av_theta_row)-log(C_beta)));
Y=(I-Mat_my-DELTA*M_XI*ones(1,N)*C_exp)\(LAM.^(-1/SIGM));
k=ones(1,N)*C_exp*Y;
k301=k;
save ('k301.mat', 'k301');


if RHO==-1
M=exp(Mat_mlamda*log(LAM)+reshape(log(GAMMA'),1,N*N)'+my*log(Y));
ks=exp(log(Y)+(I-ones(N,1)*av_theta_row')*log(LAM)+ones(N,1)*av_theta_row'*log(av_theta_row)-log(C_beta));
kcheck=(ones(1,N)*ks.^(-RHO))^(-1/RHO);
Iss=DELTA*ones(1,N)*ks;
L=exp(-log(PSSI)+log(LAM)+log(C_l)+log(Y));
X=M_XI*Iss;
nuu=exp(ln_mu)*(1-BETTA*(1-DELTA))/BETTA;
end

%Non linear System: 

if RHO<-1
    
options=optimoptions('fsolve', 'MaxFunEvals',500000, 'MaxIter', 100000);
x30=[LAM', k301];
x30=[98.5043350922051,51.6870798725288,38.8617403218898,43.2998266468549,117.118510813521,69.8761114453208,51.8500547257153,78.4017202744314,53.3500583711654,72.3292878270846,38.0539250336780,70.6091248298411,39.6562259425563,72.7614233867217,52.0540751817384,73.6791044610142,42.9893194895208,66.9769755499335,64.5753061393972,130.859869819954,41.2303863846172,39.0714457940207,43.5389230507577,28.2878832768711,26.3232395639634,32.5209589284652,16.3878344572073,29.7208148076967,16.1671817641965,10.1967886329398,0.387112206325115];
x=fsolve(@(x) n_1k_NI_eq_ss(x,BETTA,DELTA,ALFA,SIGM,GAMMA,PSSI,RHO, VARFI,1,av_theta_row),x30,options);

LAM=x(1:N)';
k=x(N+1);

CC=LAM.^(1/-SIGM);

lamda_rep=repmat(LAM,1,N)'./repmat(LAM,1,N);
Mat_my = (GAMMA.*lamda_rep); 

%Output, Materials, Sectoral Capital, Labor

Y=(I-Mat_my-DELTA*M_XI*ones(1,N)*C_exp)\(LAM.^(-1/SIGM));
k=ones(1,N)*C_exp*Y;
M=exp(Mat_mlamda*log(LAM)+reshape(log(GAMMA'),1,N*N)'+my*log(Y));
ks=exp(log(Y)+(I-ones(N,1)*av_theta_row')*log(LAM)+ones(N,1)*av_theta_row'*log(av_theta_row)-log(C_beta));
kcheck=(ones(1,N)*ks.^(-RHO))^(-1/RHO);
Iss=DELTA*ones(1,N)*ks;
L=exp(-log(PSSI)+log(LAM)+log(C_l)+log(Y));
X=M_XI*Iss;
nuu=exp(ln_mu)*(1-BETTA*(1-DELTA))/BETTA;
end

c_shrs = (LAM.^(-1/SIGM))./Y;
m_shrs = (Mat_my.*repmat(Y,1,N)')./repmat(Y,1,N);
x_shrs =X./Y;

VARFI_KSS=k^(RHO)*(VARFI.*ks'.^(-RHO));
Q_c=gma_tilde*M_c-SIGM*phi;

S_c = diag(c_shrs);
S_mm = kron(ones(1,N),m_shrs).*kron(eye(N),ones(1,N));
S_x = diag(x_shrs);

%Matrices linearized system

if SOLVE==1
nx=N+1;
ny=5*N+N^2+3;
ntot=nx+ny;
A=zeros(ntot, ntot);
B=zeros(ntot, ntot);

%percent deviations (rows for equations, col for variables)

A((N+1):(2*N),2:(N+1))=I;
A(2*N+1,6*N+N^2+3)=-BETTA*(1-DELTA);
A(2*N+1,6*N+N^2+4)=-(1-BETTA*(1-DELTA));
A(6*N+N^2+4,1)=1;


B(1:N,(N+2):(2*N+1))=S_c;
B(1:N,(3*N+2):(4*N+1))=-I;
B(1:N,(5*N+2):(5*N+N^2+1))=S_mm; 
B(1:N,(5*N+N^2+2):(6*N+N^2+1))=S_x; 

B((N+1):2*N, 2:(N+1))=VARRHO*eye(N);

B(2*N+1,6*N+N^2+3)=-1;

B(2*N+2,1)=-1;
B(2*N+2,(2*N+2):(3*N+1))=VARFI_KSS; 

B((2*N+3):(3*N+2),(2):(N+1))=I; 
B((2*N+3):(3*N+2),(2*N+2):(3*N+1))=alphad;
B((2*N+3):(3*N+2),(3*N+2):(4*N+1))=-I;
B((2*N+3):(3*N+2),(4*N+2):(5*N+1))=phi;
B((2*N+3):(3*N+2),(5*N+2):(5*N+N^2+1))=gma_tilde; 

B((3*N+3):(4*N+2),(N+2):(2*N+1))=-SIGM*I; 
B((3*N+3):(4*N+2),(3*N+2):(4*N+1))=I;
B((3*N+3):(4*N+2),(4*N+2):(5*N+1))=-I;

B((4*N+3):(4*N+N^2+2),(N+2):(2*N+1))=M_c;
B((4*N+3):(4*N+N^2+2),(3*N+2):(4*N+1))=my;
B((4*N+3):(4*N+N^2+2),(5*N+2):(5*N+N^2+1))=-eye(N*N);

B((4*N+N^2+3):(5*N+N^2+2),1)=(1+RHO)*ones(N,1);
B((4*N+N^2+3):(5*N+N^2+2),(N+2):(2*N+1))=SIGM*I; 
B((4*N+N^2+3):(5*N+N^2+2),(2*N+2):(3*N+1))=-RHO*I;
B((4*N+N^2+3):(5*N+N^2+2),(3*N+2):(4*N+1))=-I;
B((4*N+N^2+3):(5*N+N^2+2),6*N+N^2+4)=1;

B((5*N+N^2+3):(6*N+N^2+2),(N+2):(2*N+1))=SIGM*I;
B((5*N+N^2+3):(6*N+N^2+2),(5*N+N^2+2):(6*N+N^2+1))=-I;
B((5*N+N^2+3):(6*N+N^2+2),6*N+N^2+2)=ones(N,1);
B((5*N+N^2+3):(6*N+N^2+2),6*N+N^2+3)=ones(N,1);

B((6*N+N^2+3),(5*N+N^2+2):(6*N+N^2+1))=av_theta_row';
B((6*N+N^2+3),(6*N+N^2+2))=-1;

B((6*N+N^2+4),1 )=(1-DELTA);
B((6*N+N^2+4),(6*N+N^2+2))=DELTA;

%First order aproximation,The decision rules (gx) and the law of motions follow this order: 
%States: K, Zs. (law of motion wrt Kt-1, Zt-1 );Controls: Cs, ks, ys, ls, ms. (decision rule for controles wrt to states%(K, Zs)

states=N+1;
stake=1;

[gx,hx]=solab_original(A,B,states);

%x0=[k ones(N,1)'];
%[IR,IRy,IRx]=ir(gx,hx,x0,10); %Check this, it only yields the reponse of each endogenous variable to one(?) shock

end

%-----Pairwise Correlations-----%

%Data pairwise correlation output growth

data0=csvread('data_output_klems.csv');
[Tip,C]=size(data0);

data1=lagmatrix(data0,1);
yysd=(data0-data1)./data1;
dy=yysd(2:Tip,:);
dy=detrend(dy, 'constant');

[Tip,Cs]=size(dy);
dy_m=kron((mean(dy))',ones(Tip,1)');

% Covariance Matix of dy %

cov_dy=(dy-dy_m')'*(dy-dy_m')/Tip;
cov_dy_diag = diag(diag(cov_dy));
tmp=eye(N)./sqrt(cov_dy_diag);
a=diag(diag(tmp));
cor_dy=a*cov_dy*a;
cor_dy_d=cor_dy;
 for i=1:N
cor_dy(i,i)=99;
 end
 
vcor = cor_dy(triu(true(size(cor_dy))));
BB=sort(vcor,'descend');
[tot,cc]=size(BB);
vcor=BB(N+1:tot);
[fd,xid] = ksdensity(vcor);

%----One capital stock model----%

mk=hx(1,1);
M_k=mk*ones(N,1)';
m_z=hx(1,2:(N+1));

PI_ck=gx(1:N,1);
PI_cz=gx(1:N,(2:N+1));
PI_kk=gx(N+1:2*N,1);
PI_kz=gx(N+1:2*N,(2:N+1));

PI_ck=diag(PI_ck);
PI_kk=diag(PI_kk);
M_k=diag(M_k);
M_z=ones(N,1)*m_z;

B_k=PI_kk+inv(alphad)*Q_c*PI_ck;
B_z=inv(alphad)+PI_kz+inv(alphad)*Q_c*PI_cz;
    
cap_rho=B_k*M_k*inv(B_k);
Xi=B_k*(M_z-M_k*inv(B_k)*B_z);

%-------------------------------------------------------------------------%
%----------------Backing out sectoral productivity shocks------------------%
%-------------------------------------------------------------------------%

%filter using: eps: e(t)=B_z^(-1)Dy(t)-B_z^(-1)rho Dy(t-1)-B_z^(-1)Xie(t-1)

eps=zeros(Tip,N);
dy0=zeros(Tip,N);

%start with eps0=0;

dy0(2:Tip,:)=dy(2:Tip,:)*(inv(B_z))'-dy(1:Tip-1,:)*(inv(B_z)*cap_rho)';
for i=2:Tip
	eps(i,:)=dy0(i,:)-eps(i-1,:)*(inv(B_z)*Xi)';
end

eps_my=eps;

if RHO==-1  
save ('eps_my_klems_-1.mat' ,'eps_my') %saving backed out erros in case want to compare them with other models
else
save ('eps_my_klems_-2.mat' ,'eps_my')
end    

chol_my=chol(cov(eps));
save ('chol_my_klems.mat' ,'chol_my')
%chol_my=diag(diag(cov(eps)));
covmat_my=cov(eps);
save ('covmat_my_klems.mat' ,'covmat_my')

%------Carvalho's model: N capital produced with each sector output----%

%policy function +law of motion of states (From Foerster et al 2011 GAUSS code output)

M_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',1,'a1:bh60');
PII_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',2,'a1:bh90');
PIIA_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',3,'a1:ad30');

M_k=M_carvalho(1:N,1:N);
M_a=M_carvalho(1:N,N+1:2*N);

PII_ck=PII_carvalho(1:N,1:N);
PII_ca=PII_carvalho(1:N,N+1:2*N);

temp=inv(alphad);
PII_k=eye(N)+temp*Q_c*PII_ck;
PII_a=temp*(eye(N)+Q_c*PII_ca);

temp=inv(PII_k);
cap_rho=PII_k*M_k*temp;
Xi=PII_k*(M_a-M_k*temp*PII_a);

eps=zeros(Tip,N);
dy0=zeros(Tip,N);

%start with eps0=0;

dy0(2:Tip,:)=dy(2:Tip,:)*(inv(PII_a))'-dy(1:Tip-1,:)*(inv(PII_a)*cap_rho)';
for i=2:Tip
	eps(i,:)=dy0(i,:)-eps(i-1,:)*(inv(PII_a)*Xi)';
end

eps_c=eps;
save ('eps_c.mat' ,'eps_c')

covmat_c=cov(eps);
save ('covmat_c.mat' ,'covmat_c')
chol_c=chol(cov(eps));
save ('chol_c.mat' ,'chol_c')
%chol_c=diag(diag(cov(eps)));

 
%-Foerster et al 2011: N capital produced with others sectors invest. goods

M_pierre=xlsread('m_pii_30_pierre_099.xlsx',1,'a1:bh60');
PII_pierre=xlsread('m_pii_30_pierre_099.xlsx',2,'a1:bh90');
PIIA_pierre=xlsread('m_pii_30_pierre_099.xlsx',3,'a1:ad30');

M_k=M_pierre(1:N,1:N);
M_a=M_pierre(1:N,N+1:2*N);

PII_ck=PII_pierre(1:N,1:N);
PII_ca=PII_pierre(1:N,N+1:2*N);

temp=inv(alphad);
PII_k=eye(N)+temp*Q_c*PII_ck;
PII_a=temp*(eye(N)+Q_c*PII_ca);

temp=inv(PII_k);
cap_rho=PII_k*M_k*temp;
Xi=PII_k*(M_a-M_k*temp*PII_a);

eps=zeros(Tip,N);
dy0=zeros(Tip,N);

%start with eps0=0;

dy0(2:Tip,:)=dy(2:Tip,:)*(inv(PII_a))'-dy(1:Tip-1,:)*(inv(PII_a)*cap_rho)';
for i=2:Tip
	eps(i,:)=dy0(i,:)-eps(i-1,:)*(inv(PII_a)*Xi)';
end
eps_p=eps;
save ('eps_p.mat' ,'eps_p')

covmat_p=cov(eps);
save ('covmat_p.mat' ,'covmat_p')
chol_p=chol(cov(eps));
save ('chol_p.mat' ,'chol_p')
%chol_p=diag(diag(cov(eps)));

%Cholesly of variance covariance matrix of shocks

eta=vertcat( zeros(1,N),chol_my);  % cholesky for errors + zero shocks K (check this)
eta_c=vertcat( zeros(N,N), chol_c); 
eta_p=vertcat(zeros(N,N), chol_p); 

%-----------Simulated pairwise correlations----------%

 if CORRELATIONS==1

%Carvalho 2007

[A,B]=size(w');

M_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',1,'a1:bh60');
PII_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',2,'a1:bh90');
PIIA_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',3,'a1:ad30');
M_carvalho(N+1:2*N,N+1:2*N)= PII_carvalho(2*N+1:3*N,N+1:2*N); 

 %Foerster et al 2011 model

M_pierre=xlsread('m_pii_30_pierre_099.xlsx',1,'a1:bh60');
PII_pierre=xlsread('m_pii_30_pierre_099.xlsx',2,'a1:bh90');
PIIA_pierre=xlsread('m_pii_30_pierre_099.xlsx',3,'a1:ad30');
M_pierre(N+1:2*N,N+1:2*N)= PII_pierre(2*N+1:3*N,N+1:2*N); 


for i=1:SIM
e0=randn(T,N);
x0=[0 zeros(N,1)'];
x0p=[zeros(N,1)' zeros(N,1)'];

% One capital model

[YY,XX,e]=simu_1st_original(gx, hx, eta, T, x0, e0);
yys=zeros(T,N);
yys(:,:)=YY(:,(2*N+1:3*N));
kk=YY(:,(N+1:2*N));
kk1=kk(2:T,:);

%Carvalho 07

[YYc,XXc,e]=simu_1st_original(PII_carvalho(1:N,1:2*N),M_carvalho, eta_c, T, x0p, e0);
cc=YYc;
kkc=XXc(1:T,1:N); 
aa=XXc(1:T,N+1:2*N);
yysc=zeros(T,N);

for j=1:T
     yysc(j,:)=(inv(alphad)*aa(j,:)'+kkc(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end

%Foerster et al 2011

[YYp,XXp,e]=simu_1st_original(PII_pierre(1:N,1:2*N),M_pierre, eta_p, T, x0p, e0);
cc=YYp;
kkp=XXp(1:T,1:N);
aa=XXp(1:T,N+1:2*N);
yysp=zeros(T,N);

for j=1:T
yysp(j,:)=(inv(alphad)*aa(j,:)'+kkp(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end

%Look at capital now

yysl=lagmatrix(yys,1);
yys=yys-yysl;  %they are in log deviation so just y-yl es percent change 
yys=yys(2:T-1,:);

yyscl=lagmatrix(yysc,1);
yysc=yysc-yyscl;
yysc=yysc(2:T,:);

yyspl=lagmatrix(yysp,1);
yysp=yysp-yyspl;
yysp=yysp(2:T,:);

[T1,C]=size(yys);

dy=yys;
dy_m=kron((mean(yys))',ones(T1,1)');

% Covariance Matix of dy One capital model%

cov_dy=(dy-dy_m')'*(dy-dy_m')/(T1);
cov_dy_diag = diag(diag(cov_dy));
tmp=eye(N)./sqrt(cov_dy_diag);
a=diag(diag(tmp));
cor_dy=a*cov_dy*a;
cor_dy_my=cor_dy;

 for k=1:N
cor_dy(k,k)=99;
 end
 
vcor = cor_dy(triu(true(size(cor_dy))));
BB=sort(vcor,'descend');
[tot,cc]=size(BB);
vcor=BB(N+1:tot);
[fd0,xid0] = ksdensity(vcor,xid);
fmyss(:,i)=fd0;

dy=yysc;
[T1,C]=size(yysc);
dy_m=kron((mean(yysc))',ones(T1,1)');

% Covariance Matix of dy %Carvalho 07%

cov_dy=(dy-dy_m')'*(dy-dy_m')/T1;
cov_dy_diag = diag(diag(cov_dy));
tmp=eye(N)./sqrt(cov_dy_diag);
a=diag(diag(tmp));
cor_dy=a*cov_dy*a;
cor_dy_c=cor_dy;

for k=1:N
cor_dy(k,k)=99;
 end
 
vcor = cor_dy(triu(true(size(cor_dy))));
BB=sort(vcor,'descend');
[tot,cc]=size(BB);
vcor=BB(N+1:tot);
[fc0,xid0] = ksdensity(vcor,xid);
fcss(:,i)=fc0;

dy=yysp;
dy_m=kron((mean(yysp))',ones(T1,1)');

% Covariance Matix of dy Foerster et al 2011

cov_dy=(dy-dy_m')'*(dy-dy_m')/T1;
cov_dy_diag = diag(diag(cov_dy));
tmp=eye(N)./sqrt(cov_dy_diag);
a=diag(diag(tmp));
cor_dy=a*cov_dy*a;
cor_dy_p=cor_dy;

 for k=1:N
cor_dy(k,k)=99;
 end
 
vcor = cor_dy(triu(true(size(cor_dy))));
BB=sort(vcor,'descend');
[tot,cc]=size(BB);
vcor=BB(N+1:tot);
[fp0,xid0] = ksdensity(vcor,xid);
fpss(:,i)=fp0;

end

fmys=mean(fmyss,2);
fcs=mean(fcss,2);
fps=mean(fpss,2);

sector=(1:1:N);

if RHO<-1

h=figure;
plot(xid',[fd' fcs fps fmys], 'LineWidth' , 2);
title('Sectoral Output Pairwise Correlation Kernel Density')
%title('Sectoral Output Pairwise Correlation Kernel Density')
xlabel('Pairwise Correlation','LineWidth', 2) % x-axis label
ylabel('Density ','LineWidth', 2) % y-axis label
legend('Data Growth','Carvalho','Foerster et al', 'One K \rho<-1', 'One K Factor', 'Foerster et al Factor')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\sim_kernel_klems_-2.png');
%saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\sim_kernel_capital_klems_-2.png');
      
h=figure;
plot(sector',cor_dy_d(:,4),'LineWidth' , 2);
legend('Data')
hold on
plot(sector',cor_dy_c(:,4),'LineWidth' , 2);
legend('Carvalho')
hold on
plot(sector',cor_dy_p(:,4),'LineWidth' , 2);
legend('Foerster et al 2011')
hold on
plot(sector',cor_dy_my(:,4),'LineWidth' , 2);
legend('One K \rho<-1')
hold off
title('Pairwise correlations with Construction')
xlabel('sectors')
ylabel('Correlation with Construction')
legend('Data','Carvalho','Foerster et al', 'One K \rho<-1')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\corr_construction_klems_-2.png');


[fcons1,xi]=ksdensity(eps_my(:,1));
[fconsp,xi0]=ksdensity(eps_p(:,1),xi);

h=figure
plot(xi,[fconsp' fcons1'], 'LineWidth' , 2);
title('Kernel shocks to Construction')
%title('Sectoral Output Pairwise Correlation Kernel Density')
legend('FSW','One K $\rho=-1.1$')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\kernel_shock4_-2.png');


% h=figure;
% scatter(sector,diag(corr(eps_my,eps_p)),'filled', 'LineWidth', 2);
% title('Correlation backed out errors per sector One K vs Foerster et al, RHO<-1')
% xlabel('Sectors 1 to 26','LineWidth', 2) % x-axis label
% ylabel('Correlation errors different models','LineWidth', 2) % y-axis label
% saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\corr_err_p_-2.png');

end
  
if RHO==-1

h=figure;
plot(xid',[fd' fcs fps fmys], 'LineWidth' , 2);
title('Sectoral Output Pairwise Correlation Kernel Density')
%title('Sectoral Output Pairwise Correlation Kernel Density')
xlabel('Pairwise Correlation','LineWidth', 2) % x-axis label
ylabel('Density ','LineWidth', 2) % y-axis label
legend('Data Growth','Carvalho','Foerster et al', 'One K \rho=-1')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\sim_kernel_klems_-1.png');
%saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\sim_kernel_capital_klems_-1.png');
    
    
h=figure;
plot(sector',cor_dy_d(:,4),'LineWidth' , 2);
legend('Data')
hold on
plot(sector',cor_dy_c(:,4),'LineWidth' , 2);
legend('Carvalho')
hold on
plot(sector',cor_dy_p(:,4),'LineWidth' , 2);
legend('Foerster et al 2011')
hold on
plot(sector',cor_dy_my(:,4),'LineWidth' , 2);
legend('One K \rho=-1')
hold off
title('Pairwise correlations with Construction')
xlabel('sectors')
ylabel('Correlation with Construction')
legend('Data','Carvalho','Foerster et al', 'One K \rho=-1')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\corr_construction_klems_-1.png');

h=figure;
scatter(sector,diag(corr(eps_my,eps_p)),'filled', 'LineWidth', 2);
title('Correlation backed out errors per sector One K vs Foerster et al, RHO=-1')
xlabel('Sectors 1 to 26','LineWidth', 2) % x-axis label
ylabel('Correlation errors different models','LineWidth', 2) % y-axis label
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\corr_err_p_klems_-1.png');


[fcons1,xi]=ksdensity(eps_my(:,4));
[fconsp,xi0]=ksdensity(eps_p(:,4),xi);

h=figure
plot(xi,[fconsp' fcons1'], 'LineWidth' , 2);
title('Kernel shocks to Construction')
%title('Sectoral Output Pairwise Correlation Kernel Density')
legend('FSW','One K $\rho=-1$')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\kernel_shock4_-1.png');

end
end
 

%VAR Impulse Response Analysis

if VAR==1
    
%Carvalho 2007

M_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',1,'a1:bh60');
PII_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',2,'a1:bh90');
PIIA_carvalho=xlsread('m_pii_30_carvalho_099.xlsx',3,'a1:ad30');
M_carvalho(N+1:2*N,N+1:2*N)= PII_carvalho(2*N+1:3*N,N+1:2*N);

%Pierre's model
M_pierre=xlsread('m_pii_30_pierre_099.xlsx',1,'a1:bh60');
PII_pierre=xlsread('m_pii_30_pierre_099.xlsx',2,'a1:bh90');
PIIA_pierre=xlsread('m_pii_30_pierre_099.xlsx',3,'a1:ad30');
M_pierre(N+1:2*N,N+1:2*N)= PII_pierre(2*N+1:3*N,N+1:2*N);
gx_pierre=PII_pierre(1:N,1:2*N); hx_pierre=M_pierre;

BETASmy_s=zeros(N,N,SIM);
omegamy_s=zeros(N,N,SIM);

BETASc_s=zeros(N,N,SIM);
omegac_s=zeros(N,N,SIM);

BETASp_s=zeros(N,N,SIM);
omegap_s=zeros(N,N,SIM);

std_mys=zeros(N,SIM);
std_cs=zeros(N,SIM);
std_ps=zeros(N,SIM);

coms_my=zeros(N,SIM);
coms_c=zeros(N,SIM);
coms_p=zeros(N,SIM);

for i=1:SIM
e0=randn(T,N);
x0p=[zeros(N,1)' zeros(N,1)'];
x0=[0 zeros(N,1)'];

%One capital stock model

[YY,XX,e]=simu_1st_original(gx, hx, eta, T);
yys=YY(:,(2*N+1:3*N));
kk=YY(:,(N+1:2*N));
cc=YY(1:T,1:N);
aa=XX(1:T,2:N+1);

%check that this holds
for j=1:T
    yys0(j,:)=(inv(alphad)*aa(j,:)'+kk(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end

yys_a=sum(yys,2);

yysl=lagmatrix(yys,1);
yys=yys-yysl;
yys=yys(2:T,:);

yysl_a=lagmatrix(yys_a,1);
yys_a=yys_a-yysl_a;
yys_a=yys_a(2:T,:);

%Carvalho

[YYc,XXc,e]=simu_1st_original(PII_carvalho(1:N,1:2*N),M_carvalho, eta_c, T, x0p, e0);
cc=YYc;
kk=XXc(1:T,1:N);
aa=XXc(1:T,N+1:2*N);
yysc=zeros(T,N);

for j=1:T
     yysc(j,:)=(inv(alphad)*aa(j,:)'+kk(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end

yysc_a=sum(yysc,2);

%Pierre's

[YYp,XXp,e]=simu_1st_original(PII_pierre(1:N,1:2*N),M_pierre, eta_p, T, x0p, e0);
cc=YYp;
kk=XXp(1:T,1:N);
aa=XXp(1:T,N+1:2*N);
yysp=zeros(T,N);

for j=1:T
yysp(j,:)=(inv(alphad)*aa(j,:)'+kk(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end

yysp_a=sum(yysp,2);

yyscl=lagmatrix(yysc,1);
yysc=yysc-yyscl;
yysc=yysc(2:T,:);

yyscl_a=lagmatrix(yysc_a,1);
yysc_a=yysc_a-yyscl_a;
yysc_a=yysc_a(2:T,:);

yyspl=lagmatrix(yysp,1);
yysp=yysp-yyspl;
yysp=yysp(2:T,:);

yyspl_a=lagmatrix(yysp_a,1);
yysp_a=yysp_a-yyspl_a;
yysp_a=yysp_a(2:T,:);

[Y0, X] = VARmakexy(yys,1,0);
BETASmy=(X'*X)\(X'*Y0);
BETASmy_s(:,:,i)=BETASmy;

uus=Y0-X*BETASmy;
omegamy=(uus'*uus)/Tip;
omegamy_s(:,:,i)=omegamy;

[Y0, X] = VARmakexy(yysc,1,0);
BETASc=(X'*X)\(X'*Y0);
BETASc_s(:,:,i)=BETASc;

uus=Y0-X*BETASc;
omegac=(uus'*uus)/Tip;
omegac_s(:,:,i)=omegac;


[Y0, X] = VARmakexy(yysp,1,0);
BETASp=(X'*X)\(X'*Y0);
BETASp_s(:,:,i)=BETASp;

uus=Y0-X*BETASp;
omegap=(uus'*uus)/Tip;
omegap_s(:,:,i)=omegap;

std_mys(:,i)=std(yys)';
std_cs(:,i)=std(yysc)';
std_ps(:,i)=std(yysp)';

com_my=corr(yys,yys_a*ones(1,T-1));
com_c=corr(yysc,yysc_a*ones(1,T-1));
com_p=corr(yysp,yysp_a*ones(1,T-1));
coms_my(:,i)=com_my(:,1);coms_c(:,i)=com_c(:,1);coms_p(:,i)=com_p(:,1);
end

BETASmy=mean(BETASmy_s,3);
BETASp=mean(BETASp_s,3);
BETASc=mean(BETASc_s,3);

omegamy=mean(omegamy_s,3);
omegac=mean(omegac_s,3);
omegap=mean(omegap_s,3);

chol_my=chol(omegamy);
chol_c=chol(omegac);
chol_p=chol(omegap);


%coefficients cholesky transformation variance covariance matrices

chmy = omegamy(triu(true(size(omegamy))));
chmc = omegac(triu(true(size(omegac))));
chmp = omegap(triu(true(size(omegap))));

%Impulse Response Functions Models

TT=10; %time window

%IRmy=armairf({BETASmy},[],'InnovCov',chol_my,'NumObs',TT);
%IRp=armairf({BETASp},[],'InnovCov',chol_p,'NumObs',TT);
%IRc=armairf({BETASc},[],'InnovCov',chol_c,'NumObs',TT);

[IRmy] = impulse(BETASmy',omegamy,1,26,10);
[IRc] = impulse(BETASc',omegac,1,26,10);
[IRp] = impulse(BETASp',omegap,1,26,10);


%Observed VAR IRF

%Data pairwise correlation output growth

data0=csvread('data_output_klems.csv');
[Tip,C]=size(data0);
data1=lagmatrix(data0,1);
yysd=(data0-data1)./data1;
yysd=yysd(2:Tip,:);
yysd=detrend(yysd,'constant');
dy=yysd;

datak0=csvread('data_capital_klems.csv');
yysdk_a=sum(datak0,2);
dy_ak2=(yysdk_a-lagmatrix(yysdk_a,1))./lagmatrix(yysdk_a,1);
dy_ak2=dy_ak2(2:Tip,:);
dy_ak2=detrend(dy_ak2,'constant');
dy_ak=dy_ak2;

yysd_a=sum(data0,2);
dy_a2=(yysd_a-lagmatrix(yysd_a,1))./lagmatrix(yysd_a,1);
dy_a2=dy_a2(2:Tip,:);
dy_a2=detrend(dy_a2,'constant');
dy_a=dy_a2;

corr(yysdk_a, yysd_a)

%Comovements

com_my=mean(coms_my,2);
com_c=mean(com_c,2);
com_p=mean(com_p,2);
com_d=corr(dy, dy_a*ones(1,Tip));
com_d=com_d(:,1);

% Sectoral Volatilities

std_my=mean(std_mys,2);
std_c=mean(std_cs,2);
std_p=mean(std_ps,2);
std_d=std(yysd)';

[Yd, Xd] = VARmakexy(dy,1,0);
BETASd=(Xd'*Xd)\(Xd'*Yd);
uud=Yd-Xd*BETASd;
omegad=(uud'*uud)/Tip;
chol_d=chol(omegad);

[Betas,SIGMA,U,V]=olsvarc(dy,1);
chol_d=chol(SIGMA);

%Impulse Response Data

%IRd=armairf({BETASd},[], 'InnovCov',chol_d,'NumObs',TT);

%-------Check IR other codes-----%

%data=[log(g_real(st(i):Nmax1,country)) log(gdp_real(st(i):Nmax1,country)) i_sr2(st(i):Nmax1,country)];
%[N2,L2]=size(data);

%[Betas,SIGMA,U,V]=olsvarc(Yd,1);
%IRF(i,j,h) represents h-step-ahead response of variable i to 1-sigma shock of variable j. 
%[IRF, Amat, Omega] = impulse(Betas,SIGMA,1,30,10);
%[IRd] = impulse(BETASd',omegad,1,26,10);

%Weithing matrix

N2=N*N;

BETASdd=reshape(BETASd',N2,1);
BETASpp=reshape(BETASp',N2,1);
BETAScc=reshape(BETASc',N2,1);
BETASmymy=reshape(BETASmy',N2,1);

ch0 = omegad(triu(true(size(omegad))));

parameters=vertcat(BETASdd, ch0);

%-----------------Betas and Omegas-----------------%

% mmy=vertcat(BETASmymy,chmy)-vertcat(BETASdd, ch0);
% mc=vertcat(BETAScc, chmc)-vertcat(BETASdd, ch0); 
% mp=vertcat(BETASpp, chmp)-vertcat(BETASdd, ch0); 
% 
% wm=zeros(N*N+N*(N+1)/2,N*N+N*(N+1)/2);
% wm(1:N*N,1:N*N)=kron(omegad,inv(Xd'*Xd)); %Hamilton page 302 for asymp var cov of omega
% wm(N*N+1:N*N+N*(N+1)/2,N*N+1:N*N+N*(N+1)/2)=eye(N*(N+1)/2); %just identity for cholesky parameters
% 
% %GMM criteria asymp var cov omega
% 
% gmm_c=mc'*inv(wm)*mc;
% gmm_p=mp'*inv(wm)*mp;
% gmm_my=mmy'*inv(wm)*mmy;
% 
% disp('GMM criteria asymp W');
% disp([gmm_c gmm_p gmm_my]);
% 
% wm=eye(N*N+N*(N+1)/2,N*N+N*(N+1)/2);
% %wm(1:N*N,1:N*N)=kron(omegad,inv(Xd'*Xd)); %Hamilton page 302 for asymp var cov of omega
% %wm(N*N+1:N*N+N*(N+1)/2,N*N+1:N*N+N*(N+1)/2)=eye(N*(N+1)/2); %just identity for cholesky parameters
% 
% %GMM criteria identity weighting matrix
% 
% gmm_c=mc'*inv(wm)*mc;
% gmm_p=mp'*inv(wm)*mp;
% gmm_my=mmy'*inv(wm)*mmy;
% 
% disp('GMM criteria identity W');
% disp([gmm_c gmm_p gmm_my]);

%---------Only Betas---------------%

mmy=vertcat(BETASmymy)-vertcat(BETASdd);
mc=vertcat(BETAScc)-vertcat(BETASdd); 
mp=vertcat(BETASpp)-vertcat(BETASdd); 

wm=zeros(N*N,N*N);
wm(1:N*N,1:N*N)=kron(omegad,inv(Xd'*Xd)); %Hamilton page 302 for asymp var cov of omega

%GMM criteria asymp var cov omega

gmm_c=mc'*inv(wm)*mc;
gmm_p=mp'*inv(wm)*mp;
gmm_my=mmy'*inv(wm)*mmy;

disp('GMM criteria asymp W only betas');
disp([gmm_c gmm_p gmm_my]);

wm=eye(N*N,N*N);

%GMM criteria identity weighting matrix

gmm_c=mc'*inv(wm)*mc;
gmm_p=mp'*inv(wm)*mp;
gmm_my=mmy'*inv(wm)*mmy;

disp('GMM criteria identity W only betas');
disp([gmm_c gmm_p gmm_my]);

%26; power generation and gas; 22 motor; 19 equipment machinery; 4 food

t=(0:TT)';

if GRAPHS==1
    
 if RHO<-1

%Shock to Machinery

% h=figure
% plot(t,[IRd(:,4,19) IRmy(:,4,19) IRc(:,4,19) IRp(:,4,19)], 'LineWidth', 2)
% title('Response of Food to a Shock in Machinery (RHO<-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_19_4_-2.png')
% 
% h=figure
% plot(t,[IRd(:,14,19) IRmy(:,14,19) IRc(:,14,19) IRp(:,14,19)], 'LineWidth', 2)
% title('Response of Chemicals to a Shock in Machinery (RHO<-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_19_14_-2.png')


for j=1:11
ir3_19d(j) =IRd(3,19,j);     %Response of 22 to a shock in 19
ir3_19my(j) =IRmy(3,19,j);     
ir3_19c(j) =IRc(3,19,j);     
ir3_19p(j) =IRp(3,19,j);     
end


h=figure
plot(t,[ir3_19d' ir3_19my' ir3_19c' ir3_19p'], 'LineWidth', 2)
title('Response of Mining Support Activities to a Shock in Machinery (RHO=-2.1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_19_3_-2.png')

for j=1:11
ir22_19d(j) =IRd(22,19,j);     %Response of 22 to a shock in 19
ir22_19my(j) =IRmy(22,19,j);     
ir22_19c(j) =IRc(22,19,j);     
ir22_19p(j) =IRp(22,19,j);     
end

h=figure
plot(t,[ir22_19d' ir22_19my' ir22_19c' ir22_19p'], 'LineWidth', 2)
title('Response of Motor to a Shock in Machinery (RHO=-2.1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_19_22_-2.png')


% for j=1:11
% ir26_19d(j) =IRd(26,19,j);     %Response of 22 to a shock in 19
% ir26_19my(j) =IRmy(26,19,j);     
% ir26_19c(j) =IRc(26,19,j);     
% ir26_19p(j) =IRp(26,19,j);     
% end
% 
% h=figure
% plot(t,[ir26_19d' ir26_19my' ir26_19c' ir26_19p'], 'LineWidth', 2)
% title('Response of Power Generation to a Shock in Machinery (RHO<-1)')
% legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_19_26_-2.png')

% Shock to Motor
% 
% h=figure
% plot(t,[IRd(:,4,22) IRmy(:,4,22) IRc(:,4,22) IRp(:,4,22)], 'LineWidth', 2)
% title('Response of Food to a Shock in Motor (RHO<-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_22_4_-2.png')
% 
% h=figure
% plot(t,[IRd(:,14,22) IRmy(:,14,22) IRc(:,14,22) IRp(:,14,22)], 'LineWidth', 2)
% title('Response of Chemicals to a Shock in Motor (RHO<-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_22_14_-2.png')

for j=1:11
ir19_22d(j) =IRd(19,22,j);     %Response of 22 to a shock in 19
ir19_22my(j) =IRmy(19,22,j);     
ir19_22c(j) =IRc(19,22,j);     
ir19_22p(j) =IRp(19,22,j);     
end

h=figure
plot(t,[ir19_22d' ir19_22my' ir19_22c' ir19_22p'], 'LineWidth', 2)
title('Response of Machinery to a Shock in Motor (RHO=-2.1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_22_19_-2.png')

for j=1:11
ir26_22d(j) =IRd(26,22,j);     %Response of 22 to a shock in 19
ir26_22my(j) =IRmy(26,22,j);     
ir26_22c(j) =IRc(26,22,j);     
ir26_22p(j) =IRp(26,22,j);     
end

h=figure
plot(t,[ir26_22d' ir26_22my' ir26_22c' ir26_22p'], 'LineWidth', 2)
title('Response of Power Generation to a Shock in Motor (RHO=-2.1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_22_26_-2.png')

end

if RHO==-1
    
% h=figure
% plot(t,[IRd(:,4,19) IRmy(:,4,19) IRc(:,4,19) IRp(:,4,19)], 'LineWidth', 2)
% title('Response of Food to a Shock in Machinery (RHO=-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_19_4_-1.png')
% % 
% h=figure
% plot(t,[IRd(:,14,19) IRmy(:,14,19) IRc(:,14,19) IRp(:,14,19)], 'LineWidth', 2)
% title('Response of Chemicals to a Shock in Machinery (RHO=-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_19_14_-1.png')


for j=1:11
ir3_19d(j) =IRd(3,19,j);     %Response of 22 to a shock in 19
ir3_19my(j) =IRmy(3,19,j);     
ir3_19c(j) =IRc(3,19,j);     
ir3_19p(j) =IRp(3,19,j);     
end


h=figure
plot(t,[ir3_19d' ir3_19my' ir3_19c' ir3_19p'], 'LineWidth', 2)
title('Response of Mining Support Activities to a Shock in Machinery (RHO=-1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_19_3_-1.png')

for j=1:11
ir22_19d(j) =IRd(22,19,j);     %Response of 22 to a shock in 19
ir22_19my(j) =IRmy(22,19,j);     
ir22_19c(j) =IRc(22,19,j);     
ir22_19p(j) =IRp(22,19,j);     
end

h=figure
plot(t,[ir22_19d' ir22_19my' ir22_19c' ir22_19p'], 'LineWidth', 2)
title('Response of Motor to a Shock in Machinery (RHO=-1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_19_22_-1.png')


% h=figure
% plot(t,[IRd(:,26,19) IRmy(:,26,19) IRc(:,26,19) IRp(:,26,19)], 'LineWidth', 2)
% title('Response of Power Generation to a Shock in Machinery (RHO=-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_19_26_-1.png')

% Shock to Motor

% h=figure
% plot(t,[IRd(:,4,22) IRmy(:,4,22) IRc(:,4,22) IRp(:,4,22)], 'LineWidth', 2)
% title('Response of Food to a Shock in Motor (RHO=-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_22_4_-1.png')
% 
% h=figure
% plot(t,[IRd(:,14,22) IRmy(:,14,22) IRc(:,14,22) IRp(:,14,22)], 'LineWidth', 2)
% title('Response of Chemicals to a Shock in Motor (RHO=-1)')
% legend('Data','This','Carvalho','Pierre', 'location', 'northeast')
% saveas(h,'E:\Networks_in_Macroeconomics\Tex\IR_VAR_22_14_-1.png')

for j=1:11
ir19_22d(j) =IRd(19,22,j);     %Response of 22 to a shock in 19
ir19_22my(j) =IRmy(19,22,j);     
ir19_22c(j) =IRc(19,22,j);     
ir19_22p(j) =IRp(19,22,j);     
end

h=figure
plot(t,[ir19_22d' ir19_22my' ir19_22c' ir19_22p'], 'LineWidth', 2)
title('Response of Machinery to a Shock in Motor (RHO=-1)')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_22_19_-1.png')

for j=1:11
ir26_22d(j) =IRd(26,22,j);     %Response of 22 to a shock in 19
ir26_22my(j) =IRmy(26,22,j);     
ir26_22c(j) =IRc(26,22,j);     
ir26_22p(j) =IRp(26,22,j);     
end

h=figure
plot(t,[ir26_22d' ir26_22my' ir26_22c' ir26_22p'], 'LineWidth', 2)
title('Response of Power Generation to a Shock in Motor (RHO=1)')
legend('Data IP','One K model','Carvalho 2007','Foerster et al 2011', 'location', 'northeast')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\IR_VAR_22_26_-1.png')

end


%I-O graphs
labelsect={'Oil&gas','Min','Min supp.','Food','Tobacco','Text Mills','Text Prod','Apparel','Leather','Wood','Pulp','Printing','Petrol.','Chemicals','Plastic','Non Met. Min.','Metal. Min.','Fab. Metals','Machinery','Computer','Elect. Equip.','Motor','Furniture','Medical Equip.','Logging','Power Gen.'};

io_usa=(log(GAMMA+0.0000001));
fpath = 'C:\Users\uqjmira1\Dropbox\Miranda-Young\Comparing_Models\';

h=figure;
imagesc(sqrt(GAMMA));
colorbar;
title('U.S. Industrial Sectors I-O 1997')
xlabel('Sectors Buying Input')
ylabel('Sectors Supplying Input')
set(gca,'YLim',[1 26]);
set(gca,'YTick',[1:26]);
set(gca,'YTickLabel',labelsect)
h1 = colorbar;
title(h1,'Sqrt Int. Input Share')
saveas(h, fullfile(fpath, 'usaiocolor97'), 'png')


load theta97
theta_usa=(log(THETA+0.0000001));

h=figure;
imagesc(sqrt(THETA));
colorbar;
title('U.S. Ind. Sectors Capital Flow 1997')
xlabel('Sectors Buying Invest. goods')
ylabel('Sectors Supplying Invest. goods')
set(gca,'YLim',[1 26]);
set(gca,'YTick',1:26);
set(gca,'YTickLabel',labelsect)
h1 = colorbar;
title(h1,'Sqrt Invest. Goods Share')
saveas(h, fullfile(fpath, 'usacapflowcolor97'), 'png')
toc

end
end
