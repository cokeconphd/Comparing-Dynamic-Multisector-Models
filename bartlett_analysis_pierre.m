%Multisector intratemporal adjustment costs of capital 

%clearvars
clc  
tic

% Key values

shry = xlsread('1997DetailedItemOutput.xlsx',4,'f3:f28' );
sig_ee = xlsread('sig_ee_factor_pierre.xlsx',1,'a1:z26' ); %variance covariance shocks in Pierre's structural factor model
sig_chol_ee=chol(sig_ee); 

T=1000; %size of the sample
SIM=30; %number of simulation to then average
w=0:0.01:2.14;
nimp=1000;

% Parameter Values

BETTA = 0.99;      %discount rate
DELTA = 0.025;     %depreciation rate
SIGM  = 1;         %CRRA coeficient Utility function
RHO   = -1;        %Intratemporal adjustment cost of capital parameter (-1 means no cost)
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

GAMMA=GAMA;
clear GAMA %To avoid problems of names 

for i=1:col
    for j=1:col
        if GAMMA(i,j)==0
    GAMMA(i,j) = 10e-20;
        end
    end
end

% Parameter/Matrix Values

Mat_mlamda = kron(ones(N,1), eye(N)) - kron(eye(N),ones(N,1)); 
my=kron(ones(N,1),eye(N));
alphad = diag(ALFA,0);
phi = eye(N) - alphad - repmat(sum(GAMMA)', 1, N).*eye(N);
gma_tilde0 = (kron(GAMMA',ones(1,N)));
gma_tilde1 = kron(ones(1,N),eye(N));
gma_tilde  =  gma_tilde0.*gma_tilde1;
z=ones(N,1);
C_l = phi*ones(N,1);
C_beta = ((1-BETTA*(1-DELTA))./(BETTA*ALFA)); 
I=eye(N);
BETTA_TILDE=1-BETTA*(1-DELTA);

M_c=SIGM*(kron(eye(N),ones(N,1)) - kron(ones(N,1),eye(N))); 
 
Q_c=gma_tilde*M_c-SIGM*phi;


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

load chol_c
load chol_p
load covmat_p
load covmat_c

eta_c=vertcat(chol_c, zeros(N,N));
eta_p=vertcat(chol_p, zeros(N,N));

%-----------Spectral Analysis------------%


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

[A,B]=size(w');

s_y=zeros(N,N,A);


for j=1:A
A0=(eye(N)-cap_rho*exp(-1i*w(j)))*(eye(N)-cap_rho*exp(1i*w(j))).';
B0=(eye(N)+Xi*exp(-1i*w(j)))*(PII_a)'*covmat_c*PII_a*(eye(N)+Xi*exp(1i*w(j))).';
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
B0=(eye(N)+Xi*exp(-1i*w(j)))*(PII_a)'*covmat_p*PII_a*(eye(N)+Xi*exp(1i*w(j))).';
s_y(:,:,j)=(2*pi)^(-1)*(B0*inv(A0));
end

varcov=sum(s_y,3);

s_y_agg_p2=zeros(A,1);

for j=1:A
s_y_agg_p2(j)=shry'*(s_y(:,:,j)./varcov)*shry;
end

% Barlett Estimates

s_y_agg_c_s=zeros(A,SIM);
s_y_agg_p_s=zeros(A,SIM);

for i=1:SIM
e0=randn(T,N);
x0p=[zeros(N,1)' zeros(N,1)'];

%Carvalho

[YYc,XXc,e]=simu_1st_original(PII_carvalho(1:N,1:2*N),M_carvalho, eta_c, T, x0p, e0);
cc=YYc;
kk=XXc(1:T,1:N);
aa=XXc(1:T,N+1:2*N);
yysc=zeros(T,N);

for j=1:T
     yysc(j,:)=(inv(alphad)*aa(j,:)'+kk(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end


%Pierre's

[YYp,XXp,e]=simu_1st_original(PII_pierre(1:N,1:2*N),M_pierre, eta_p, T, x0p, e0);
cc=YYp;
kk=XXp(1:T,1:N);
aa=XXp(1:T,N+1:2*N);
yysp=zeros(T,N);

for j=1:T
yysp(j,:)=(inv(alphad)*aa(j,:)'+kk(j,:)'+inv(alphad)*Q_c*cc(j,:)')';
end

% yyscl=lagmatrix(yysc,1);
% yysc=yysc-yyscl; %or
% yysc=(yysc-yyscl).*y;
% yysc=yysc(2:T,:);

% yyspl=lagmatrix(yysp,1);
% yysp=yysp-yyspl;
% yysp=yysp(2:T,:);

%Carvalho
[T1,aver]=size(yysp);

yhat = mean(yysc); % mean of columns of matrix X
yhat =kron(yhat,ones(1,T1)');
Yc0=yysc-yhat; % centered data matrix
Sigma0 = (1/T1)*(Yc0'*Yc0); % covariance matrix

q=20;

sigmas=zeros(N,N,q);
Ycs=zeros(T1,N,q);

for l=1:q
    Ycs(:,:,l)=lagmatrix(yysc,l)-yhat;
    sigmas(:,:,l)=(1/T1)*Yc0(l+1:T1,:)'*Ycs(l+1:T1,:,l);
end

s_y_c=zeros(N,N,A);

for j=1:A
  for k=1:q
    if k==1
    s_y_c(:,:,j)=Sigma0;
    end
  s_y_c(:,:,j)=s_y_c(:,:,j)+ (1-k/(q+1))*(sigmas(:,:,k)*exp(-1i*k*w(j))+sigmas(:,:,k)'*exp(1i*k*w(j)));
  end
end

varcov=sum(s_y_c,3);

s_y_agg_c=zeros(A,1);

for j=1:A
s_y_agg_c(j)=shry'*(s_y_c(:,:,j)./varcov)*shry;
end

%Pierre

yhat = mean(yysp); % mean of columns of matrix X
yhat =kron(yhat,ones(1,T1)');
Yc0=yysp-yhat; % centered data matrix
Sigma0 = (1/T1)*(Yc0'*Yc0); % covariance matrix

q=20;

sigmas=zeros(N,N,q);
Ycs=zeros(T1,N,q);

for l=1:q
    Ycs(:,:,l)=lagmatrix(yysp,l)-yhat;
    sigmas(:,:,l)=(1/T1)*Yc0(l+1:T1,:)'*Ycs(l+1:T1,:,l);
end

s_y_p=zeros(N,N,A);

for j=1:A
  for k=1:q
    if k==1
    s_y_p(:,:,j)=Sigma0;
    end
  s_y_p(:,:,j)=s_y_p(:,:,j)+ (1-k/(q+1))*(sigmas(:,:,k)*exp(-1i*k*w(j))+sigmas(:,:,k)'*exp(1i*k*w(j)));
  end
end

varcov=sum(s_y_p,3);

s_y_agg_p=zeros(A,1);

for j=1:A
s_y_agg_p(j)=shry'*(s_y_p(:,:,j)./varcov)*shry;
end

s_y_agg_c_s(:,i)=s_y_agg_c;
s_y_agg_p_s(:,i)=s_y_agg_p;

end

s_y_aggc=mean(s_y_agg_c_s,2);
s_y_aggp=mean(s_y_agg_p_s,2);

% toc

h=figure;
plot(w,[s_y_agg_c2 s_y_aggc],'LineWidth', 2);
legend('Actual', 'Bartlett q=20','location','northeast')
title('Carvalho 2007 Spectrum Aggregate Output Growth')
%saveas(h,'E:\Networks_in_Macroeconomics\Tex\bartlett_c.png');

h=figure;
plot(w,[s_y_agg_p2 s_y_aggp],'LineWidth', 2);
legend(' Actual','Bartlett q=20','location','northeast')
title('Foerester et al 2011 Spectrum Aggregate Output Growth')
%saveas(h,'E:\Networks_in_Macroeconomics\Tex\bartlett_p.png');


h=figure;
plot(w,[s_y_agg_c2 s_y_agg_p2],'LineWidth', 2);
legend(' Actual','Bartlett q=20','location','northeast')
title('Foerester et al 2011 Spectrum Aggregate Output Growth')
%saveas(h,'E:\Networks_in_Macroeconomics\Tex\bartlett_p.png');



