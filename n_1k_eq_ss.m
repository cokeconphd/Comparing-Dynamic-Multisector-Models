function y = n_1k_eq_ss(x,BETTA,DELTA,ALFA,SIGM,GAMMA,PSSI,RRHO,VARFI,nk)

%Defining variables

[row,col]=size(GAMMA);

LAM=zeros(row,1)';

for i=1:row
LAM(i)=x(i);
end

%mus
LAMN=ones(row,1)*LAM(nk);

%k bar
K=x(row+1);
Kbar=zeros(row,1);
Kbar(nk)=K;
MKbar=zeros(row,1);
MKbar(nk)=1;
KKK=ones(row,1)*K;

%II=DELTA*K;
%IIN=[zeros(row-1,1)' II]';

Mat_mlamda = kron(ones(row,1), eye(row)) - kron(eye(row),ones(row,1)); 
my=kron(ones(row,1),eye(row));
alphad = diag(ALFA,0);
phi = eye(row) - alphad - repmat(sum(GAMMA)', 1, row).*eye(row);
gma_tilde0 = (kron(GAMMA',ones(1,row)));
gma_tilde1 = kron(ones(1,row),eye(row));
gma_tilde  =  gma_tilde0.*gma_tilde1;
z=ones(row,1);

%------------System of Equations--------------%

ialrho=inv(alphad*(RRHO+1)/RRHO);
C_l = phi*ones(row,1);
C_beta = ((1-BETTA*(1-DELTA))./(BETTA*ALFA)); 
lamda_rep=repmat(LAM',1,row)'./repmat(LAM',1,row);
Mat_my = (GAMMA.*lamda_rep); 
I=eye(row);
B=I+alphad/RRHO-repmat(sum(GAMMA)', 1, row).*eye(row)-phi;


%N EQUATIONS (EQ 39 APPENDIX)

LLAMMDAS=log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*MKbar*K))-inv(B)*z-inv(B)*alphad./RRHO*(log(LAMN)-log(LAM')+log(VARFI')+(1+RRHO)*ones(row,1)*log(K)+log(C_beta))-inv(B)*(gma_tilde*(Mat_mlamda*log(LAM')+reshape(log(GAMMA'),1,row*row)') + phi*(-log(PSSI)+log(LAM')+log(C_l)));
%LLAMMDAS=log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+ELTA*MKbar*K))-ialrho*z-ialrho*alphad./RRHO*(log(LAMN)-log(LAM')+log(VARFI')+(1+RRHO)*ones(row,1)*log(K)+log(C_beta))-ialrho*(gma_tilde*(Mat_mlamda*log(LAM')+reshape(log(GAMMA'),1,row*row)') + phi*(-log(PSSI)+log(LAM')+log(C_l)));

for i=1:row
y(i) = LLAMMDAS(i);
end


%1 EQUATION (EQUATION 40 APPENDIX FROM MKT CLEARIN OF CAPITAL)
y(row+1)=K-(VARFI*(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar).*(LAMN.^(-1).*C_beta.^(-1).*LAM'.*KKK.^(-1-RRHO))))^(-1/RRHO);

%y(row+1)=K-(exp((-1/RRHO)*(log(VARFI)-log(LAMN)+log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))-log(C_beta)+log(LAM')-(1+RRHO)*log(KKK))));
%y(row+1)=K+(exp((1/RRHO)*(log(VARFI)-log(LAMN)+log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))-log(C_beta)+log(LAM')-(1+RRHO)*log(KKK))))^(-1/RRHO);
%y(row+1)=K-(VARFI*(exp((1/RRHO)*(log(LAMN)-log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))+log(C_beta)-log(LAM')+log(VARFI')+(1+RRHO)*log(KKK)))).^(-RRHO))^(-1/RRHO);
%y(row+1)=K-([0.1, 0.3, 0.4, 0.1, 0.5]*(exp((1/RRHO)*(log(LAMN)-log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))+log(C_beta)-log(LAM')+(1+RRHO)*log(KKK)))).^(-RRHO))^(-1/RRHO);

