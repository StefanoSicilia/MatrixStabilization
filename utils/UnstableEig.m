function [lambda,Leig,Reig,gamma]=UnstableEig(M,m,delta)
%% UnstableEig:
% Computes the eigenvalues of M with real part larger than -d<0, where
% d=delta or d=delta(2). It is assumed that these eigenvalues are p<=m.
% The function also returns the matrices Leig and Reig with the left and 
% right unit eigenvectors x_i and y_i as columns associated to \lambda_i.
% Finally the vector gamma contains the real scalars
% gamma_i=(real(\lambda_i(M))+\delta)/(x_i'*y_i) if length(delta)=1
% gamma_i=rhod_i*(dphi_i*rhod_i+2*phi_i)/(x_i'*y_i) if length(delta)=2,
% where rhod=real(lambda+delta(2)) and
% [phi,dphi]=Psi_HI(delta(1),d,rhod),
% with x_i and y_i normalized such that x_i'*y_i>0.
    
    dl=length(delta);
    switch dl
        case 1
            d=delta;
        case 2
            d=delta(2);
        otherwise
            error('Incorrect length for delta')
    end
    
    %% MAIN COMPUTATION
    n=size(M,1);
    m=min(m,n);
    [Reig,Rlambda]=eigs(M,m,'largestreal');
    [Leig,Llambda]=eigs(M.',m,'largestreal');  
    [lambda,ind1]=sort(diag(Rlambda),'descend','ComparisonMethod','real');
    [~,ind2]=sort(diag(Llambda),'descend','ComparisonMethod','real');
    Reig=Reig(:,ind1);
    Leig=conj(Leig(:,ind2));
    j=m;
    real_lambda=real(lambda);
    while (j>1) && (real_lambda(j)<-d)
        j=j-1;
    end
    mstar=j;
    
    %% TARGET EIGENVALUES AND EIGENVECTORS FOUND
    lambda=lambda(1:mstar);
    Reig=Reig(:,1:mstar);
    Leig=Leig(:,1:mstar);
        
    %% COMPUTATION OF SCALAR PRODUCTS
    scalars=zeros(mstar,1);
    for i=1:mstar
       scalars(i)=Leig(:,i)'*Reig(:,i); 
    end
    
    %% ROTATION OF THE EIGENVECTORS
    absscalars=abs(scalars);
    Reig=Reig./(ones(n,1)*(absscalars./scalars)');
    
    %% EVENTUAL CONDITION NUMBER
    if nargout>3
        switch length(delta)
            case 1
                gamma=real(lambda+d)./absscalars;
            case 2
                rho=real(lambda);
                [phi,dphi]=Psi_HI(-d,-delta(1),rho);
                gamma=(rho+d).*(dphi.*(rho+d)+2*phi)./(2*absscalars);
        end
    end
end