%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-function that computes impulse response (not confidence intervals) 
function [IRF, Amat, Omega] = impulse(Betas,Omega,nlag,ndim,hmax)
          % Omega: covariance matrix of error term of (p,1)-autoregression 
          % Cholesky decomposition. P is a lower triangular such that Omega = PP'.
          
          Omega=Omega(1:ndim,1:ndim);
          P = chol(Omega, 'lower');
          Betas=Betas(1:ndim,:)';
          % # of lags: nlag = result_VAR.nlag;
          % dimension of VAR: ndim = result_VAR.ndim;

          % construct a ndim x ndim x nlag matrix of OLS estimates from result_VAR.OLSE1.
          % Note how result_VAR.OLSE1 is constructed.
          % It is a (nlag * ndim) x ndim matrix containing least squares estimates for
          % (nlag,1)-autoregression.
          % The first ndim x ndim block is the TRANSPOSE of A_1 and so on.
          % Amat removes the transpose.
          Amat = zeros(ndim, ndim, nlag);
          for j = 1:nlag
               first = (j-1)*ndim + 1;
               last = j*ndim;
               Amat(:,:,j) = Betas(first:last, :)';
          end;
           % Construct VMA(infty) coefficients up to lag hmax
          % Recursive formula: Psi(j) = sum_{s=1}^{nlag} A(s) Psi(j-s)         
          Psimat = zeros(ndim, ndim, hmax);
          indmat = repmat((1:hmax)', 1, nlag) - repmat(1:nlag, hmax, 1);  % index matrix
          % indmat = [   1 - 1, ... ,     1 - nlag;
          %                :  ,  :  ,       :     ;
          %           hmax - 1, ... ,  hmax - nlag];
          
          for j = 1:hmax   % compute Psi(j)
               for l = 1:nlag      % add nlag terms
                    ind = indmat(j,l);
                    if ind < 0     % Psi with negative index is a null matrix
                        increment = zeros(ndim, ndim);
                    elseif ind == 0   % Psi with zero index is eye(K)
                        increment = Amat(:,:,l);
                    else           % Psi with positive index is just as is
                        increment = Amat(:,:,l) * Psimat(:,:,ind);
                    end;    
                    Psimat(:,:,j) = Psimat(:,:,j) + increment;
               end;     
          end;    

          % compute impulse response functions
          % IRF(i,j,i) represents i-step-ahead response of variable i to 
          % 1-sigma shock of variable j. 
          IRF = zeros(ndim, ndim, hmax+1);
          IRF(:,:,1) = P;  % 0-step ahead IRF
          for i = 1:hmax
               IRF(:, :, i+1) = Psimat(:,:,i) * P;
          end;    
end  %   close sub-