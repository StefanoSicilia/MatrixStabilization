function [U,S,V,F,Times,Ranks,derfeps]=InnerIter_SA(A,epsilon,opts)
%% InnerIter_SA: Structured Inner Iteration with adaptive-rank integrator
% Performs the inner iteration of the two-level procedure for determining
% the closest stable matrix to A implemented by Stabilize.
% STRUCTURED distance with ADAPTIVE-RANK integrator.

    %% Inizialization   
    % Variables
    n=size(A,1);
    fun=opts.fun;
    delta=opts.delta;
    maxit=opts.maxit;
    h0=opts.h;
    tol=opts.tol;
    theta=opts.theta;
    safestop=opts.safestop;
    r=opts.r; 
    ranktol=opts.ranktol;
    row=opts.row;
    col=opts.col;
    
    % Function and output structures
    switch fun
        case 'standard'
            feval=@(x) 0.5*(sum((x+delta).^2));
        case 'Hermite'
            feval=@(x) 0.5*sum(((x+delta).^2).*Psi_HI(-2*delta,-delta,x));
            delta=[delta,2*delta];
    end
    F=zeros(maxit,1);
    Times=zeros(maxit,1);
    Ranks=zeros(maxit,1);
    
    % Eventual starting point (if not provided in the input)
    if ~isfield(opts,'initialguess')
        [lambda,Leig,Reig,gamma]=UnstableEig(A,r,delta);
        f=feval(real(lambda));
        [U,TL]=qr(Leig,0);
        [V,TR]=qr(Reig,0);
        S=-TL*diag(gamma)*TR';
    else
        InitialGuess=opts.initialguess;
        U=InitialGuess.U0;
        S=InitialGuess.S0;
        V=InitialGuess.V0;     
        f=-1;
        gamma=[];
    end 
    E=projsparse(row,col,U,S,V);
    normE=norm(E,'fro');
    E=E/normE;
    S=S/normE;
    F(1)=f;
    Ranks(1)=length(gamma);

    %% First iteration
    [lambda,Leig,Reig,gamma]=UnstableEig(A+epsilon*E,r,delta);
    Gamma=diag(gamma);
    GS=projsparse(row,col,Leig,Gamma,Reig);    
    mu=real(trace((U'*(GS*V))*S));
    nu=norm(GS,'fro');
    f=feval(real(lambda));
    g=epsilon*(nu^2-mu^2);
    hj=h0; 
    F(2)=f;
    Times(2)=hj;
    Ranks(2)=length(gamma);
    j=2;

    %% Main computations
    while j<=maxit && g>tol && Ranks(j)>0
        fh=f;
        cont=1;
        h=hj;
        while fh>=f && cont<safestop
            
            % Initial and augmented rank 
            r0=size(S,1);
            rho=min([2*r0,n]);

            % Update U,S and V
            eta=ComputEta(Leig,Gamma,Reig,U,V,E);
            alpha=1+h*eta;
    
            % i) K-step
            K0=U*S;
            K1=alpha*K0-h*GS*V;
            [Uhat,~]=qr([K1,U],0);
            Mhat=Uhat'*U;

            % ii) L-step
            L0=V*S';
            L1=alpha*L0-h*(GS')*U;
            [Vhat,~]=qr([L1,V],0);
            Nhat=Vhat'*V;

            % iii) S-step
            Shat0=Mhat*S*Nhat';
            Ehat0=projsparse(row,col,Uhat,Shat0,Vhat);
            normEhat0=norm(Ehat0,'fro');
            Ehat0=Ehat0/normEhat0;
            Shat0=Shat0/normEhat0;
            [~,Leig,Reig,gamma]=UnstableEig(A+epsilon*Ehat0,rho,delta);
            Gamma=diag(gamma);
            etahat0=ComputEta(Leig,Gamma,Reig,Uhat,Vhat,Ehat0);
            Shat1=(1+h*etahat0)*Shat0-h*(Uhat'*Leig)*Gamma*(Reig'*Vhat);
            
            % iv) Truncation
            [Phat,Sigmahat,Qhat]=svd(Shat1);
            sumsig=0;
            s=diag(Sigmahat).^2;
            tolsq=ranktol^2;
            r1=rho;
            while sumsig<tolsq 
                sumsig=sumsig+s(r1);
                r1=r1-1;
            end
            P1=Phat(:,1:r1);
            Q1=Qhat(:,1:r1);

            % v) Fixed-rank integrator results
            S1=Sigmahat(1:r1,1:r1);
            Uh=Uhat*P1;
            Vh=Vhat*Q1;
            Eh=projsparse(row,col,Uh,S1,Vh);
            normEh=norm(Eh,'fro');
            Eh=Eh/normEh;
            Sh=S1/normEh;
            
            % New eigenvalues and gradient
            [lambdah,Leigh,Reigh,gammah]=UnstableEig(A+epsilon*Eh,r1,delta);
            Gammah=diag(gammah);
            GSh=projsparse(row,col,Leigh,Gammah,Reigh);
            muh=real(trace((Uh'*(GSh*Vh))*Sh));
            nuh=norm(GSh,'fro');
            fh=feval(real(lambdah));
            gh=epsilon*(nuh^2-muh^2);

            % Eventual enlargement of the step
            if fh>=f
                h=h/theta;
            end
            cont=cont+1;
        end

        % Armijo Stepsize Selection
        if cont==safestop 
            g=tol;
            j=j+1;
        else
            hnext=h;
            if fh>=f-(h/theta)*g
               hnext=h/theta; 
            end

            if hnext==hj %%

                % Initial and augmented rank
                r0=size(S,1);
                rho=min([2*r0,n]);
                
                % Update U,S and V with stepsize ht=h*theta
                ht=h*theta;                        
                alpha=1+ht*eta;

                % i) K-step
                K0=U*S;
                K1=alpha*K0-ht*GS*V;
                [Uhat,~]=qr([K1,U],0);
                Mhat=Uhat'*U;

                % ii) L-step
                L0=V*S';
                L1=alpha*L0-ht*(GS')*U;
                [Vhat,~]=qr([L1,V],0);
                Nhat=Vhat'*V;

                % iii) S-step
                Shat0=Mhat*S*Nhat';
                Ehat0=projsparse(row,col,Uhat,Shat0,Vhat);
                normEhat0=norm(Ehat0,'fro');
                Ehat0=Ehat0/normEhat0;
                Shat0=Shat0/normEhat0;
                [~,Leig,Reig,gamma]=UnstableEig(A+epsilon*Ehat0,rho,delta);
                Gamma=diag(gamma);
                eta=ComputEta(Leig,Gamma,Reig,Uhat,Vhat,Ehat0);
                Shat1=(1+ht*eta)*Shat0-ht*(Uhat'*Leig)*Gamma*(Reig'*Vhat);
                
                % iv) Truncation
                [Phat,Sigmahat,Qhat]=svd(Shat1);
                sumsig=0;
                s=diag(Sigmahat).^2;
                tolsq=ranktol^2;
                r1=rho;
                while sumsig<tolsq 
                    sumsig=sumsig+s(r1);
                    r1=r1-1;
                end
                P1=Phat(:,1:r1);
                Q1=Qhat(:,1:r1);

                % v) Fixed-rank integrator results
                S1=Sigmahat(1:r1,1:r1);
                Ut=Uhat*P1;
                Vt=Vhat*Q1;
                Et=projsparse(row,col,Ut,S1,Vt);
                normEt=norm(Et,'fro');
                Et=Et/normEt;
                St=S1/normEt;

                % New eigenvalues and gradient
                [lambdat,Leigt,Reigt,gammat]=UnstableEig(A+epsilon*Et,r1,delta);
                Gammat=diag(gammat);
                GSt=Leigt*Gammat*Reigt';
                mut=real(trace((Ut'*(GSt*Vt))*St));
                nut=norm(GSt,'fro');
                ft=feval(real(lambdat));
                gt=epsilon*(nut^2-mut^2);

                if fh>ft
                    Uh=Ut;
                    Sh=St;
                    Vh=Vt;
                    Eh=Et;
                    Leigh=Leigt;
                    Gammah=Gammat;
                    Reigh=Reigt;
                    GSh=GSt;
                    nuh=nut;
                    fh=ft;
                    gh=gt;
                    hnext=ht;
                end
            end

            % Update for the next iteration
            U=Uh;
            S=Sh;
            V=Vh;
            E=Eh;
            Leig=Leigh;
            Gamma=Gammah;
            Reig=Reigh;
            GS=GSh;
            nu=nuh;
            f=fh;
            g=gh;
            hj=hnext; 
            j=j+1;
            F(j)=f;
            Times(j)=hj;
            Ranks(j)=size(S,1);
        end
    end 
    derfeps=-nu;

    %% FINAL OUTPUTS
    iterations=max(j-1,2);
    F=F(1:iterations);
    Times=cumsum(Times(1:iterations));
    Ranks=Ranks(1:iterations);

end