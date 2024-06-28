function [U,S,V,F,Times,Ranks,derfeps]=InnerIter_UF(A,epsilon,opts)
%% InnerIter_UF: Unstructured Inner Iteration with fixed rank integration
% Performs the inner iteration of the two-level procedure for determining
% the closest stable matrix to A implemented by Stabilize.
% UNSTRUCTURED distance with FIXED-RANK integrator.

    %% Inizialization   
    % Variables
    fun=opts.fun;
    delta=opts.delta;
    maxit=opts.maxit;
    h0=opts.h;
    tol=opts.tol;
    theta=opts.theta;
    safestop=opts.safestop;
    r=opts.r; 
    
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
    E=U*S*V';
    normE=norm(E,'fro');
    E=E/normE;
    S=S/normE;
    F(1)=f;
    Ranks(1)=length(gamma);

    %% First iteration
    [lambda,Leig,Reig,gamma]=UnstableEig(A+epsilon*E,r,delta);
    G=Leig*diag(gamma)*Reig';    
    mu=real(trace((U'*(G*V))*S));
    nu=norm(G,'fro');
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

            % Update U,S and V            
            alpha=1+h*mu;
    
            % i) K-step
            K0=U*S;
            K1=alpha*K0-h*G*V;
            [Uhat,~]=qr(K1,0);
            Mhat=Uhat'*U;

            % ii) L-step
            L0=V*S';
            L1=alpha*L0-h*(G')*U;
            [Vhat,~]=qr(L1,0);
            Nhat=Vhat'*V;

            % iii) S-step
            Shat0=Mhat*S*Nhat';
            Ehat0=Uhat*Shat0*Vhat';
            [~,Leig,Reig,gamma]=UnstableEig(A+epsilon*Ehat0,r,delta);
            Gamma=diag(gamma);
            Ghat0=Leig*Gamma*Reig';
            mu=real(trace((Vhat'*Reig)*(Gamma*(Leig'*Uhat)*Shat0)));
            Shat1=(1+mu*h)*Shat0-h*Uhat'*Ghat0*Vhat;

            % iv) Fixed-rank integrator results
            [P1,S1,Q1]=svd(Shat1);
            Sh=S1/norm(S1,'fro');
            Uh=Uhat*P1;
            Vh=Vhat*Q1;
            
            % New eigenvalues and gradient
            Eh=Uh*Sh*Vh';
            [lambdah,Leigh,Reigh,gammah]=UnstableEig(A+epsilon*Eh,r,delta);
            Gammah=diag(gammah);
            Gh=Leigh*Gammah*Reigh';
            muh=real(trace((Uh'*(Gh*Vh))*Sh));
            nuh=norm(Gh,'fro');
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

                % Update U,S and V with stepsize ht=h*theta
                ht=h*theta;                        
                alpha=1+ht*mu;

                % i) K-step
                K0=U*S;
                K1=alpha*K0-ht*G*V;
                [Uhat,~]=qr(K1,0);
                Mhat=Uhat'*U;

                % ii) L-step
                L0=V*S';
                L1=alpha*L0-ht*(G')*U;
                [Vhat,~]=qr(L1,0);
                Nhat=Vhat'*V;

                % iii) S-step
                Shat0=Mhat*S*Nhat';
                Ehat0=Uhat*Shat0*Vhat';
                [~,Leig,Reig,gamma]=UnstableEig(A+epsilon*Ehat0,r,delta);
                Gamma=diag(gamma);
                Ghat0=Leig*Gamma*Reig';
                mu=real(trace((Vhat'*Reig)*(Gamma*(Leig'*Uhat)*Shat0)));
                Shat1=(1+mu*ht)*Shat0-ht*Uhat'*Ghat0*Vhat;

                % iv) Fixed-rank integrator results
                [P1,S1,Q1]=svd(Shat1);
                St=S1/norm(S1,'fro');
                Ut=Uhat*P1;
                Vt=Vhat*Q1;

                % New eigenvalues and gradient
                Et=Ut*St*Vt';
                [lambdat,Leigt,Reigt,gammat]=UnstableEig(A+epsilon*Et,r,delta);
                Gammat=diag(gammat);
                Gt=Leigt*Gammat*Reigt';
                mut=real(trace((Ut'*(Gt*Vt))*St));
                nut=norm(Gt,'fro');
                ft=feval(real(lambdat));
                gt=epsilon*(nut^2-mut^2);

                if fh>ft
                    Uh=Ut;
                    Sh=St;
                    Vh=Vt;
                    nuh=nut;
                    fh=ft;
                    gh=gt;
                    hnext=ht;
                    muh=mut;
                end
            end

            % Update for the next iteration
            U=Uh;
            S=Sh;
            V=Vh;
            nu=nuh;
            f=fh;
            g=gh;
            hj=hnext; 
            mu=muh;
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