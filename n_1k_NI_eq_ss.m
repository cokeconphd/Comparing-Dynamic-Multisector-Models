function y = n_1k_NI_eq_ss(x,BETTA,DELTA,ALFA,SIGM,GAMMA,PSSI,RRHO,VARFI,nk, theta)

%Defining variables

[N,col]=size(GAMMA);

LAM=zeros(N,1)';

for i=1:N
LAM(i)=x(i);
end

K=x(N+1);
KKK=ones(N,1)*K;

Mat_mlamda = kron(ones(N,1), eye(N)) - kron(eye(N),ones(N,1)); 
my=kron(ones(N,1),eye(N));
alphad = diag(ALFA,0);
phi = eye(N) - alphad - repmat(sum(GAMMA)', 1, N).*eye(N);
gma_tilde0 = (kron(GAMMA',ones(1,N)));
gma_tilde1 = kron(ones(1,N),eye(N));
gma_tilde  =  gma_tilde0.*gma_tilde1;
z=ones(N,1);
ln_mu=theta'*log(LAM')-theta'*log(theta);
M_XI=exp(ln_mu)*theta./LAM';

%------------System of Equations--------------%

C_l = phi*ones(N,1);
C_beta = ((1-BETTA*(1-DELTA))./(BETTA*ALFA)); 
lamda_rep=repmat(LAM',1,N)'./repmat(LAM',1,N);
Mat_my = (GAMMA.*lamda_rep); 
I=eye(N);
B=I+alphad/RRHO-repmat(sum(GAMMA)', 1, N).*eye(N)-phi;


LLAMMDAS=log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*M_XI*K))-inv(B)*z-inv(B)*alphad./RRHO*((ones(N,1)*theta'-I)*log(LAM')-ones(N,1)*theta'*log(theta)+(1+RRHO)*ones(N,1)*log(K)+log(C_beta))-inv(B)*(gma_tilde*(Mat_mlamda*log(LAM')+reshape(log(GAMMA'),1,N*N)') + phi*(-log(PSSI)+log(LAM')+log(C_l)));

for i=1:N
y(i) = LLAMMDAS(i);
end

y(N+1)=K-(ones(1,N)*(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*M_XI*K).*C_beta.^(-1).*KKK.^(-1-RRHO).*(exp((1/RRHO)*((ones(N,1)*theta'-I)*log(LAM')-ones(N,1)*theta'*log(theta)))).^(-RRHO)))^(-1/RRHO);

%y(N+1)=K-(exp((-1/RRHO)*(log(VARFI)-log(LAMN)+log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))-log(C_beta)+log(LAM')-(1+RRHO)*log(KKK))));
%y(N+1)=K+(exp((1/RRHO)*(log(VARFI)-log(LAMN)+log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))-log(C_beta)+log(LAM')-(1+RRHO)*log(KKK))))^(-1/RRHO);
%y(N+1)=K-(VARFI*(exp((1/RRHO)*(log(LAMN)-log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))+log(C_beta)-log(LAM')+log(VARFI')+(1+RRHO)*log(KKK)))).^(-RRHO))^(-1/RRHO);
%y(N+1)=K-([0.1, 0.3, 0.4, 0.1, 0.5]*(exp((1/RRHO)*(log(LAMN)-log(inv(I-Mat_my)*(LAM'.^(-1/SIGM)+DELTA*Kbar))+log(C_beta)-log(LAM')+(1+RRHO)*log(KKK)))).^(-RRHO))^(-1/RRHO);

