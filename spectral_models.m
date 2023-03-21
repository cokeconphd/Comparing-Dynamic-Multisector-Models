%Multisector intratemporal adjustment costs of capital 

clearvars
clc  
tic

% Key values

shry = xlsread('1997DetailedItemOutput.xlsx',4,'f3:f28' );
sig_ee = xlsread('sig_ee_factor_pierre.xlsx',1,'a1:z26' ); %variance covariance shocks in Pierre's structural factor model
sig_chol_ee=chol(sig_ee); 

T=500; %size of the sample
SIM=30; %number of simulation to then average
w=0:0.01:2.14;
nimp=1000;
figures=1;

% Parameter Values

BETTA = 0.99;      %discount rate
DELTA = 0.025;     %depreciation rate
SIGM  = 1;         %CRRA coeficient Utility function
RHO   = -2.1;        %Intratemporal adjustment cost of capital parameter (-1 means no cost)
VARRHO= 1;         %persistence o productivity shock
SIGZ  =0.0076;    %volatility of technology shock (this will have to be a matriz later on
PSSI  = 1;        %Coeficient for labor/leiusure utility

%Load the IO shares and production function shares

load gamma26
load alfa26
load av_sell_theta26

%capital good sector

NK=19;

[N, col]=size(GAMA);

%adjusting theta vector

rm=0;
trick=zeros(N,1);
trick(19,:)=rm;
av_theta1 = av_sell_theta26 + trick;
av_theta1 = av_theta1./(1+rm);
av_theta_row=av_theta1;

GAMMA=GAMA;
clear GAMA %To avoid problems of names 

for i=1:col
    for j=1:col
        if GAMMA(i,j)==0
    GAMMA(i,j) = 10e-20;
        end
    end
end

% Steady State and Parameter/Matrix Values

VARFI=ones(1,N);
Mat_mlamda = kron(ones(N,1), eye(N)) - kron(eye(N),ones(N,1)); 
my=kron(ones(N,1),eye(N));
alphad = diag(ALFA,0);
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
LAM261=LAM;
save ('LAM261.mat','LAM261');
CC=LAM.^(1/-SIGM);

lamda_rep=repmat(LAM,1,N)'./repmat(LAM,1,N);
Mat_my = (GAMMA.*lamda_rep); 

C_exp=diag(exp((I-ones(N,1)*av_theta_row')*log(LAM)+ones(N,1)*av_theta_row'*log(av_theta_row)-log(C_beta)));
Y=(I-Mat_my-DELTA*M_XI*ones(1,N)*C_exp)\(LAM.^(-1/SIGM));
k=ones(1,N)*C_exp*Y;
k261=k;
save ('k261.mat', 'k261');


if RHO==-1
M=exp(Mat_mlamda*log(LAM)+reshape(log(GAMMA'),1,N*N)'+my*log(Y));
ks=exp(log(Y)+(I-ones(N,1)*av_theta_row')*log(LAM)+ones(N,1)*av_theta_row'*log(av_theta_row)-log(C_beta));
kcheck=(ones(1,N)*ks.^(-RHO))^(-1/RHO);
Iss=DELTA*ones(1,N)*ks;
L=exp(-log(PSSI)+log(LAM)+log(C_l)+log(Y));
X=M_XI*Iss;
nuu=exp(ln_mu)*(1-BETTA*(1-DELTA))/BETTA;
end

% %Non linear System: Initial Values

if RHO<-1
    
load LAM261
load k261
%x26=[LAM261' k261];
options=optimoptions('fsolve', 'MaxFunEvals',500000, 'MaxIter', 100000);
x26=[1.65063895821368,2.37133813148524,2.03844708736734,1.94138942740135,1.62530407204538,3.58118381674838,2.87741596326143,2.56356418466516,3.37163712683570,3.13611867025373,3.66012344810571,2.47584239656941,1.67536043294499,2.48983902702571,2.78072820461841,2.06622243021480,4.00387804299062,2.47044020951162,4.16767612567084,2.02065508209771,3.61510992291013,6.42185998790056,3.13654048009798,2.52951020971712,1.00773675348846,1.12211465385685,18.3844903180997];
x=fsolve(@(x) n_1k_NI_eq_ss(x,BETTA,DELTA,ALFA,SIGM,GAMMA,PSSI,RHO, VARFI,NK,av_theta_row),x26,options);

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

M_clev=-SIGM*diag(M)*(Mat_mlamda/diag(CC));
M_ylev=diag(M)*(my*diag(Y.^(-1)));
M_mlev=eye(N^2);
gma_tilde_lev =(gma_tilde/diag(M)).*repmat(diag(Y), 1, N); %check
c_shrs = (LAM.^(-1/SIGM))./Y;
Q_c=gma_tilde*M_c-SIGM*phi;
m_shrs0 = (Mat_my.*repmat(Y,1,N)')./repmat(Y,1,N);
m_shrs = kron(ones(1,N),m_shrs0).*kron(eye(N),ones(1,N));
VARFI_KSS=k^(1+RHO)*(VARFI.*ks'.^(-RHO-1));

shry26=sum(Y);
shry26=Y/shry26;

%Matrices Log linearized system

nx=N+1;
ny=5*N+N^2+3;
ntot=nx+ny;
A=zeros(ntot, ntot);
B=zeros(ntot, ntot);

%level deviations (rows for equations, col for variables)

A((N+1):(2*N),2:(N+1))=I;
A(2*N+1,6*N+N^2+3)=-BETTA*(1-DELTA);
A(2*N+1,6*N+N^2+4)=-BETTA;
A(6*N+N^2+4,1)=1;


B(1:N,(N+2):(2*N+1))=I;
B(1:N,(3*N+2):(4*N+1))=-I;
B(1:N,(5*N+2):(5*N+N^2+1))=S_m; 
B(1:N,(5*N+N^2+2):(6*N+N^2+1))=I; 

B((N+1):2*N, 2:(N+1))=VARRHO*eye(N);

B(2*N+1,6*N+N^2+3)=-1;

B(2*N+2,1)=-1;
B(2*N+2,(2*N+2):(3*N+1))=VARFI_KSS; %ok

B((2*N+3):(3*N+2),(2):(N+1))=diag(Y).*I*(exp(1)^(-1)); %ok
B((2*N+3):(3*N+2),(2*N+2):(3*N+1))=alphad.*diag(Y).*diag(ks.^(-1));%ok
B((2*N+3):(3*N+2),(3*N+2):(4*N+1))=-I;%ok
B((2*N+3):(3*N+2),(4*N+2):(5*N+1))=phi.*diag(Y).*diag(L.^(-1));%ok
B((2*N+3):(3*N+2),(5*N+2):(5*N+N^2+1))=gma_tilde_lev; %ok

B((3*N+3):(4*N+2),(N+2):(2*N+1))=-SIGM*I.*diag(CC.^(-1-SIGM))*phi*diag(Y); %ok
B((3*N+3):(4*N+2),(3*N+2):(4*N+1))=I*diag(CC.^(-SIGM))*phi;
B((3*N+3):(4*N+2),(4*N+2):(5*N+1))=-I*PSSI;

B((4*N+3):(4*N+N^2+2),(N+2):(2*N+1))=M_clev;
B((4*N+3):(4*N+N^2+2),(3*N+2):(4*N+1))=M_ylev;
B((4*N+3):(4*N+N^2+2),(5*N+2):(5*N+N^2+1))=-M_mlev;

B((4*N+N^2+3):(5*N+N^2+2),1)=(1+RHO)*ones(N,1)*nuu*k^(RHO);
B((4*N+N^2+3):(5*N+N^2+2),(N+2):(2*N+1))=SIGM*I.*diag(CC.^(-1 - SIGM)).*alphad.*diag(Y).*diag(ks.^(RHO)); %or nuu*k^(1+RHO)*(SIGM)*diag(CC.^(-1))
B((4*N+N^2+3):(5*N+N^2+2),(2*N+2):(3*N+1))=-RHO*nuu*k^(1+RHO)*diag(ks.^(-1)); %or -RHO*nuu*k^(1+RHO)*diag(ks.^(-1))
B((4*N+N^2+3):(5*N+N^2+2),(3*N+2):(4*N+1))=-nuu*k^(1+RHO)*diag(Y.^(-1));
B((4*N+N^2+3):(5*N+N^2+2),6*N+N^2+4)=k^(1+RHO);

B((5*N+N^2+3):(6*N+N^2+2),(N+2):(2*N+1))=SIGM*I.*diag(CC.^(-1-SIGM));
B((5*N+N^2+3):(6*N+N^2+2),(5*N+N^2+2):(6*N+N^2+1))=-diag(av_theta_row*exp(ln_mu)./(X.^2));
B((5*N+N^2+3):(6*N+N^2+2),6*N+N^2+2)=av_theta_row*exp(ln_mu)./X;
B((5*N+N^2+3):(6*N+N^2+2),6*N+N^2+3)=av_theta_row*Iss./X;

B((6*N+N^2+3),(5*N+N^2+2):(6*N+N^2+1))=(av_theta_row*Iss./X)';
B((6*N+N^2+3),(6*N+N^2+2))=-1;

B((6*N+N^2+4),1 )=(1-DELTA);
B((6*N+N^2+4),(6*N+N^2+2))=1;

%First order aproximation,The decision rules (gx) and the law of motions follow this order: 
%States: K, Zs. (law of motion wrt Kt-1, Zt-1 );Controls: Cs, ks, ys, ls, ms. (decision rule for controles wrt to states%(K, Zs)

states=N+1;
stake=1;

[gx,hx]=solab_original(A,B,states)



%-----Pairwise Correlations-----%

data0=xlsread('ipl2.xlsx',2,'c2:ab145');
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

 for i=1:N
cor_dy(i,i)=99;
 end
 
vcor = cor_dy(triu(true(size(cor_dy))));
BB=sort(vcor,'descend');
[tot,cc]=size(BB);
vcor=BB(N+1:tot);
[fd,xid] = ksdensity(vcor);

%Intratemp Adj cost model

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

Q_c_my=gma_tilde_lev*M_clev-SIGM*phi*diag(Y./CC);
B_k=PI_kk*diag(Y./ks)+inv(alphad)*Q_c_my*PI_ck;
B_z=inv(alphad)*diag(Y./exp(1))+PI_kz*diag(Y./ks)+inv(alphad)*Q_c_my*PI_cz;
    
cap_rho=B_k*M_k*inv(B_k);
Xi=B_k*(M_z-M_k*inv(B_k)*B_z);

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
save ('eps_my_-1.mat' ,'eps_my')
else
save ('eps_my_-2.mat' ,'eps_my')
end    

chol_my=chol(cov(eps));
save ('chol_my.mat' ,'chol_my')

covmat_my=cov(eps);
save ('covmat_my.mat' ,'covmat_my')
se=std(eps)'; %Get the covariance matrix no just the diagonal
se_mat=diag(se.^2);
cov_dx = zeros(N,N);
   
 i=1;
 while (i) ;  
    x=B_z(:,i)*se(i);
    cov_dx=cov_dx+(x*x');
 t=1;
 while (t);
 x=cap_rho*x;
    if t==1;
        x=x+Xi(:,i)*se(i);
    end;
 cov_dx=cov_dx+(x*x');
 t=t+1;
     if t>nimp
        break;
     end;
 end
 i=i+1;
    if i>N
    break
    end
 end

cov_dx_diag = diag(diag(cov_dx));
tmp=eye(N)./sqrt(cov_dx_diag);
a=diag(diag(tmp));
cor_dx=a*cov_dx*a;

 for i=1:N
cor_dx(i,i)=99;
 end
 
vcorx = cor_dx(triu(true(size(cor_dx))));
BB=sort(vcorx,'descend');
[tot,cc]=size(BB);
vcorx=BB(N+1:tot);
[fmy,ximy] = ksdensity(vcorx,xid);


%-----------Spectral Analysis------------%

load chol_my
load chol_p
load covmat_my
[A,B]=size(w');

%Exact spectral density (varma(1,1)
    
s_y=zeros(N,N,A);

omega=diag(diag(covmat_my)); %only sectoral shocks

for j=1:A
A0=(eye(N)-cap_rho*exp(-1i*w(j)))*(eye(N)-cap_rho*exp(1i*w(j))).';
B0=(eye(N)+Xi*exp(-1i*w(j)))*(B_z)'*omega*B_z*(eye(N)+Xi*exp(1i*w(j))).';
s_y(:,:,j)=(2*pi)^(-1)*(B0*inv(A0));
end

varcov=sum(s_y,3);

s_y_agg=zeros(A,1);

for j=1:A
s_y_agg(j)=shry'*(s_y(:,:,j)./varcov)*shry;
end


%Carvalho 2007

M_carvalho=xlsread('m_pii_26_carvalho.xlsx',1,'a1:az52');
PII_carvalho=xlsread('m_pii_26_carvalho.xlsx',2,'a1:az78');
PIIA_carvalho=xlsread('m_pii_26_carvalho.xlsx',3,'a1:z26');
M_carvalho(N+1:2*N,N+1:2*N)= PII_carvalho(2*N+1:3*N,N+1:2*N);

% %Pierre's model
M_pierre=xlsread('m_pii_26_pierre.xlsx',1,'a1:az52');
PII_pierre=xlsread('m_pii_26_pierre.xlsx',2,'a1:az78');
PIIA_pierre=xlsread('m_pii_26_pierre.xlsx',3,'a1:z26');
M_pierre(N+1:2*N,N+1:2*N)= PII_pierre(2*N+1:3*N,N+1:2*N);


[A,B]=size(w');
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

s_y=zeros(N,N,A);

for j=1:A
A0=(eye(N)-cap_rho*exp(-1i*w(j)))*(eye(N)-cap_rho*exp(1i*w(j))).';
B0=(eye(N)+Xi*exp(-1i*w(j)))*(PII_a)'*omega*PII_a*(eye(N)+Xi*exp(1i*w(j))).';
s_y(:,:,j)=(2*pi)^(-1)*(B0*inv(A0));
end

varcov=sum(s_y,3);

s_y_agg_c2=zeros(A,1);

for j=1:A
s_y_agg_c2(j)=shry'*(s_y(:,:,j)./varcov)*shry;
end

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

[A,B]=size(w');

s_y=zeros(N,N,A);

for j=1:A
A0=(eye(N)-cap_rho*exp(-1i*w(j)))*(eye(N)-cap_rho*exp(1i*w(j))).';
B0=(eye(N)+Xi*exp(-1i*w(j)))*(PII_a)'*omega*PII_a*(eye(N)+Xi*exp(1i*w(j))).';
s_y(:,:,j)=(2*pi)^(-1)*(B0*inv(A0));
end

varcov=sum(s_y,3);

s_y_agg_p2=zeros(A,1);

for j=1:A
s_y_agg_p2(j)=shry'*(s_y(:,:,j)./varcov)*shry;
end

save ('s_y_agg_p2.mat','s_y_agg_p2');
save ('s_y_agg_c2.mat','s_y_agg_c2');


if RHO==-1
s_y_agg1=s_y_agg;
save ('s_y_agg.mat','s_y_agg');
end

if RHO==-1.5
s_y_agg15=s_y_agg;
save ('s_y_agg15.mat','s_y_agg15');
end

if RHO==-2.1
s_y_agg17=s_y_agg;
save ('s_y_agg17.mat','s_y_agg17');
end


kk=100;

if figures==1
load s_y_agg17
load s_y_agg15
load s_y_agg    
h=figure;
plot(w(1:kk),[s_y_agg(1:kk) s_y_agg15(1:kk) s_y_agg17(1:kk) s_y_agg_c2(1:kk) s_y_agg_p2(1:kk)],'LineWidth', 2);
legend('One K \rho=-1','One K \rho=-1.5','One K \rho=-2.1', ' Carvalho', 'Foerster et al','location','northeast')
title('Spectrum Aggregate Output Growth')
xlabel('Frequency \omega')
saveas(h,'C:\Users\uqjmira1\Dropbox\Miranda-Young\Matching_Sectoral_Dynamics\spectral_all.png');
end

