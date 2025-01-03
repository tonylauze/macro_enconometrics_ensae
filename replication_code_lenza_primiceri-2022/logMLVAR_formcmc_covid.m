function [logML,betadraw,drawSIGMA]=logMLVAR_formcmc_covid(par,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,draw,hyperpriors,priorcoef,Tcovid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the log-posterior (or the logML if hyperpriors=0),
% and draws from the posterior distribution of the coefficients and of the 
% covariance matrix of the residuals of the BVAR of Giannone, Lenza and 
% Primiceri (2015), augmented with a change in volatility at the time of 
% Covid (March 2020).
%
% Last modified: 06/02/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=par(1);

d=n+2;

theta=MIN.theta;
miu=MIN.miu;
eta=MIN.eta';

psi=SS*(d-n-1);

if isempty(Tcovid);
    if sur==1;
        theta=par(2);
        if noc==1;
            miu=par(3);
        end
    elseif sur==0;
        if noc==1;
            miu=par(2);
        end
    end
else
    ncp=4;      % number of covid parameters
    eta=par(2:ncp+1);
    
    invweights = ones(T,1); 
    invweights(Tcovid)=eta(1);
    invweights(Tcovid+1)=eta(2);
    if T>Tcovid+1;
        invweights(Tcovid+2:T)=1+(eta(3)-1)*eta(4).^[0:T-Tcovid-2];
    end
    y = diag(1./invweights)*y;
    x = diag(1./invweights)*x;
    
    if sur==1;
        theta=par(ncp+2);
        if noc==1;
            miu=par(ncp+3);
        end
    elseif sur==0;
        if noc==1;
            miu=par(ncp+2);
        end
    end
end

if mn.alpha==0;
    alpha=2;
elseif mn.alpha==1;
    alpha=par(end);
end


%% return a very low value of the posterior if the parameters are outside the bounds
if sum([lambda;eta;theta;miu;alpha]<[MIN.lambda;MIN.eta';MIN.theta;MIN.miu;MIN.alpha])>0 | ...
        sum([lambda;eta;theta;miu]>[MAX.lambda;MAX.eta';MAX.theta;MAX.miu])>0
    logML=-10e15;
    betadraw=[];
    drawSIGMA=[];
    return
else
    
    %% priors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=1+n*lags;
    omega=zeros(k,1);
    omega(1)=Vc;
    for i=1:lags
        omega(1+(i-1)*n+1:1+i*n)=(d-n-1)*(lambda^2)*(1/(i^alpha))./psi;
    end
    
    % prior scale matrix for the covariance of the shocks
    PSI=diag(psi);
    
    Td=0;
    xdsur=[];
    ydsur=[];
    xdnoc=[];
    ydnoc=[];
    % dummy observations if sur and/or noc = 1
    if sur==1;
        
        xdsur=[1/theta (1/theta)*repmat(y0,1,lags)];
        ydsur=(1/theta)*y0;
        
        y=[y;ydsur];
        x=[x;xdsur];
        Td=1;
    end
    
    if noc==1;
        
        ydnoc=(1/miu)*diag(y0); ydnoc(pos,pos)=0;
        xdnoc=[zeros(n,1) (1/miu)*repmat(diag(y0),1,lags)];
        
        y=[y;ydnoc];
        x=[x;xdnoc];
        Td=Td+n;
    end
    
    
    %% output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % posterior mode of the VAR coefficients
    betahat=(x'*x+diag(1./omega))\(x'*y+diag(1./omega)*b);
    
    % VAR residuals
    epshat=y-x*betahat;
    
    % logML
    T=T+Td;
    
    aaa=diag(sqrt(omega))*(x'*x)*diag(sqrt(omega));
    bbb=diag(1./sqrt(psi))*(epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b))*diag(1./sqrt(psi));
    
    eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
    eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
    
    logML = - n*T*log(pi)/2 + sum(gammaln((T+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) +...
        - T*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
    
    if sur==1 | noc==1;
        yd=[ydsur;ydnoc];
        xd=[xdsur;xdnoc];
        
        % prior mode of the VAR coefficients
        % betahatd=(xd'*xd+diag(1./omega))\(xd'*yd+diag(1./omega)*b);
        betahatd=b;     % this is the case for our priors (the line above delivers the same but is numerically not very stable)
        
        
        % VAR residuals at the prior mode
        epshatd=yd-xd*betahatd;
        
        aaa=diag(sqrt(omega))*(xd'*xd)*diag(sqrt(omega));
        bbb=diag(1./sqrt(psi))*(epshatd'*epshatd + (betahatd-b)'*diag(1./omega)*(betahatd-b))*diag(1./sqrt(psi));
        
        eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
        eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
        
        % normalizing constant
        norm = - n*Td*log(pi)/2 + sum(gammaln((Td+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) +...
            - Td*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
        
        logML=logML-norm;
    end
    
    % log determinant to account for the re-weighting
    if ~isempty(Tcovid)
        logML = logML - n*sum(log(invweights));
    end
    
    if hyperpriors==1;
        logML=logML+logGammapdf(lambda,priorcoef.lambda.k,priorcoef.lambda.theta);
        if sur==1;
            logML=logML+logGammapdf(theta,priorcoef.theta.k,priorcoef.theta.theta);
        end
        if noc==1;
            logML=logML+logGammapdf(miu,priorcoef.miu.k,priorcoef.miu.theta);
        end
        if ~isempty(Tcovid)            
            logML = logML - 2*log(eta(1)) - 2*log(eta(2)) - 2*log(eta(3)) +...
                    logBetapdf(eta(4),priorcoef.eta4.alpha,priorcoef.eta4.beta);
        end
    end
    
    % takes a draw from the posterior of SIGMA and beta, if draw is on
    if draw==0
        betadraw=[];
        drawSIGMA=[];
    elseif draw==1
        S=PSI + epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b);
        
        [V E]=eig(S);
        Sinv=V*diag(1./abs(diag(E)))*V';
        eta=mvnrnd(zeros(1,n),Sinv,T+d);
        drawSIGMA=(eta'*eta)\eye(n);
        %[cholSIGMA,junk]=chol(drawSIGMA);
        %betadraw=betahat+mvnrnd(zeros(k,1),(x'*x+diag(1./omega))\eye(k),n)'*cholSIGMA;
        cholSIGMA=cholred((drawSIGMA+drawSIGMA')/2);
        cholZZinv = cholred((x'*x+diag(1./omega))\eye(k));
        betadraw=betahat + cholZZinv'*randn(size(betahat))*cholSIGMA;
    end
end

function r=logGammapdf(x,k,theta);
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logBetapdf(x,al,bet);
r=(al-1)*log(x)+(bet-1)*log(1-x)-betaln(al,bet);

function r=logIG2pdf(x,alpha,beta);
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);

function C = cholred(S);
[v,d] = eig((S+S')/2);
d = diag(real(d));
warning off
scale = mean(diag(S))*1e-12;
J = (d>scale);
C = zeros(size(S));
C(J,:) = (v(:,J)*(diag(d(J)))^(1/2))';