function [d,U,S,V,FinalF,outiter]=Stabilize(A,InnerPar,OuterPar)
%% Stabilize: 
% Given an unstable matrix A, it computes an approximation of the 
% (un)structured distance to stability and the associated perturbation.
% InnerPar and OuterPar are structs containing the parameters for the inner
% and outer iterations respectively and OuterPar.type specifies of distance
% sought and the integrator to use:
%   -'UA': Unstructured and Adaptive rank,
%   -'UF': Unstructured and Fixed rank,
%   -'SA': Structured and Adaptive rank,
%   -'SF': Structured and Fixed rank.
% The distance found is d, while the perturbed and stable matrix A+d*E is
% described by the unit Forbenius norm matrix E, where 
%   -E=U*S*V' for the unstructured distance,
%   -E=proj(U*S*V') for the structured distance, where proj projects onto
%   the set of real matrices with the same sparsity path of A.
%
% Inputs
%   - A: an unstable sparse matrix,
%   - InnerPar: inner iteration parameters stored in a struct with fields
%       - delta, ensures strict stability
%       - r, upper bound of the unstable eigenvalues of A
%       - fun, describes the objective function ('standard' or 'Hermite')
%       - maxit, maximum number of inner iterations
%       - tol, tolerance on the error
%       - h, initial stepsize for the integrator
%       - theta, parameter for the Armijo rule
%       - safestop, parameter to ensure monotonicity of the functional
%       - ranktol, tolerance for the rank in the adaptive integrator
%       - row, rows of the pattern of A
%       - col, columns of the pattern of A.
%   - OuterPar: outer iteration parameters stored in a struct with fields
%       - epsilon=[el,eu,e0], lower bound, upper bound and initial guess
%           for the distance
%       - niter, maximum number of outer iterations
%       - toler, tolerance on the error.
%
% For some examples of choice of the parameters see the numerical
% experiments in the folder 'tests'.
%
% It implements the two-level method proposed in
% N. Guglielmi, S.Sicilia, "Stabilization of a matrix via a 
% low-rank-adaptive ODE", https://arxiv.org/html/2402.14657v1

    %% Inizialization       
    % Inner variables
    delta=InnerPar.delta;
    r=InnerPar.r;
    fun=InnerPar.fun;
    switch fun
        case 'standard'
            feval=@(x) 0.5*(sum((x+delta).^2));
        case 'Hermite'
            feval=@(x) 0.5*sum(((x+delta).^2).*Psi_HI(-2*delta,-delta,x));
    end
    
    % Outer variables
    el=OuterPar.epsilon(1);
    eu=OuterPar.epsilon(2);
    e0=OuterPar.epsilon(3);
    niter=OuterPar.niter;
    toler=OuterPar.toler;    
    
    %% Definition of the method and the starting point
    switch OuterPar.type
        case 'UF'
            % Choice of the Inner Iteration
            InnerIter=@(epsilon,InnerPar) InnerIter_UF(A,epsilon,InnerPar);
            
            % Starting point
            [~,Leig,Reig,gamma]=UnstableEig(A,r,delta);
            [rowLeig,colLeig]=size(Leig);
            Leig=[Leig,zeros(rowLeig,r-colLeig)];
            Reig=[Reig,zeros(rowLeig,r-colLeig)];
            gamma=[gamma;zeros(r-colLeig,1)];
            [U,TL]=qr(Leig,0);
            [V,TR]=qr(Reig,0);
            S=-TL*diag(gamma)*TR';
            E=U*S*V';
            normE=norm(E,'fro');
            E=E/normE;
            S=S/normE;
            
            % Starting gradient
            [lambda,Leig,Reig,gamma]=UnstableEig(A+e0*E,r,delta);
            f=feval(real(lambda));
            GS=Leig*diag(gamma)*Reig';
            fprime=-norm(GS,'fro');
        case 'SF'
            % Choice of the Inner Iteration
            InnerIter=@(epsilon,InnerPar) InnerIter_SF(A,epsilon,InnerPar);
            row=InnerPar.row;
            col=InnerPar.col;
            
            % Starting point
            [~,Leig,Reig,gamma]=UnstableEig(A,r,delta);
            [rowLeig,colLeig]=size(Leig);
            Leig=[Leig,zeros(rowLeig,r-colLeig)];
            Reig=[Reig,zeros(rowLeig,r-colLeig)];
            gamma=[gamma;zeros(r-colLeig,1)];
            [U,TL]=qr(Leig,0);
            [V,TR]=qr(Reig,0);
            S=-TL*diag(gamma)*TR';
            E=projsparse(row,col,U,S,V);
            normE=norm(E,'fro');
            E=E/normE;
            S=S/normE;
            
            % Starting gradient
            [lambda,Leig,Reig,gamma]=UnstableEig(A+e0*E,r,delta);
            f=feval(real(lambda));
            GS=projsparse(row,col,Leig,diag(gamma),Reig);
            fprime=-norm(GS,'fro');
        case 'UA'
             % Choice of the Inner Iteration
            InnerIter=@(epsilon,InnerPar) InnerIter_UA(A,epsilon,InnerPar);
            
            % Starting point
            [~,Leig,Reig,gamma]=UnstableEig(A,r,delta);
            [U,TL]=qr(Leig,0);
            [V,TR]=qr(Reig,0);
            S=-TL*diag(gamma)*TR';
            E=U*S*V';
            normE=norm(E,'fro');
            E=E/normE;
            S=S/normE;
            
            % Starting gradient
            [lambda,Leig,Reig,gamma]=UnstableEig(A+e0*E,r,delta);
            f=feval(real(lambda));
            GS=Leig*diag(gamma)*Reig';
            fprime=-norm(GS,'fro');
        case 'SA'
            % Choice of the Inner Iteration
            InnerIter=@(epsilon,InnerPar) InnerIter_SA(A,epsilon,InnerPar);
            row=InnerPar.row;
            col=InnerPar.col;
            
            % Starting point
            [~,Leig,Reig,gamma]=UnstableEig(A,r,delta);
            [U,TL]=qr(Leig,0);
            [V,TR]=qr(Reig,0);
            S=-TL*diag(gamma)*TR';
            E=projsparse(row,col,U,S,V);
            normE=norm(E,'fro');
            E=E/normE;
            S=S/normE;
            
            % Starting gradient
            [lambda,Leig,Reig,gamma]=UnstableEig(A+e0*E,r,delta);
            f=feval(real(lambda));
            GS=projsparse(row,col,Leig,diag(gamma),Reig);
            fprime=-norm(GS,'fro');
    end    
    
    %% Newton-Bisection method
    j=1;
    epsilonvec=zeros(niter,1);
    epsilonvec(1)=e0;
    epsilon=e0;
    while j<=niter && (eu-el)>toler
        if f<toler
            eu=min(eu,epsilon);
            epsilon=0.5*(eu+el);
        else
            el=max(el,epsilon);
            epsilon=epsilon-f/fprime;
        end
        if epsilon<el || epsilon>eu
            epsilon=0.5*(el+eu);
        end
        InnerPar.initialguess.U0=U;
        InnerPar.initialguess.S0=S;
        InnerPar.initialguess.V0=V;
        [U,S,V,F,~,~,fprime]=InnerIter(epsilon,InnerPar);
        f=F(end);
        j=j+1;
        epsilonvec(j)=epsilon;
    end

    %% Final Output
    outiter=j-1;
    d=epsilonvec(outiter);
    FinalF=f;
    disp(['Number of outer iterations: ',num2str(outiter),...
        ' out of ',num2str(niter),'.'])



end