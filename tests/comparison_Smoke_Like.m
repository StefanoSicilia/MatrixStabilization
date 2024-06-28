%% Comparison between Stabilize, Gillis-Sharma and Noferini-Poloni methods
% The script requires the following algorithms:
% 1) the method developed by Gillis and Sharma in  
% N. Gillis and P. Sharma. "On computing the distance to stability for 
% matrices using linear dissipative Hamiltonian systems",
% and implemented by the function 'neareststablealgo', whose code can be
% found in
% https://sites.google.com/site/nicolasgillis/code
% 2) the method developed by Noferini and Poloni in
% V. Noferini and F. Poloni. "Nearest \Omega-stable matrix via Riemannian 
% optimization", 
% and implemented by the function 'nearest_stable_complex' from the 
% repository
% https://github.com/fph/nearest-omega-stable
    
    %% Example generation
    n=30;
    B=gallery('smoke',n);
    rng(1)
    [Q,~]=qr(randn(n)+1i*randn(n));
    A=Q*B*Q';
    r=15;
    v=sort(eig(B),'descend','ComparisonMethod','real');
    A=sparse(A);
    
    %% Parameters of the inner iteration
    delta=0;
    maxit=150;
    h=0.1; 
    tol=1e-9;
    theta=1.3;
    safestop=15;
    stdfun='standard';
    ranktol=1e-9;
    
    %% Parameters of the outer iteration
    el=1e-3;
    eu=10;
    e0=1e-1;
    epsilon=[el,eu,e0];
    niter=200;
    toler=tol;
    
    %% Grouping of the parameters
    InnerParsStd=struct('delta',delta,'maxit',maxit,'h',h,'tol',tol,...
        'theta',theta,'r',r,'safestop',safestop,'fun',stdfun,...
        'ranktol',ranktol);
    OuterPar=struct('epsilon',epsilon,'niter',niter,'toler',toler);
    
    %% Structures for the outputs and plots 
    % Sizes of samples
    nex=1;
    totalsize=nex+4;
    
    % Structures for all results
    Dist=zeros(totalsize,1);
    Ranks=zeros(totalsize,1);
    Fun=zeros(totalsize,1);
    Abscissa=zeros(totalsize,1);
    CPUtimes=zeros(totalsize,1);
    eig_tot={[2,1]};
    eig_stable={[2,1]};
    eig_unclear={[2,1]};
    eig_unstable={[2,1]};
    
    %% Main Computation: Standard functional
    disp('Standard------------------------------')
    
    % Adaptive
    OuterPar.type='UA';
    j=1;
    tic;
    [Dist(j),U,S,V,Fun(j)]=Stabilize(A,InnerParsStd,OuterPar);
    CPUtimes(j)=toc;
    eig_tot{j}=eig(full(A+Dist(j)*U*S*V'));
    eig_stable{j}=eig_tot{j}(eig_tot{j}<-delta);
    eig_unclear{j}=eig_tot{j}(eig_tot{j}<=0 & eig_tot{j}>=-delta);
    eig_unstable{j}=eig_tot{j}(eig_tot{j}>0);
    Abscissa(j)=max(real(eig_tot{j}));
    Ranks(j)=size(S,1);
    
    % Fixed rank
    OuterPar.type='UF';
    for j=2:nex
        rk=r+j-2;
        InnerParsStd.r=rk;
        tic;
        [Dist(j),U,S,V,Fun(j)]=Stabilize(A,InnerParsStd,OuterPar);
        CPUtimes(j)=toc;
        Abscissa(j)=max(real(eig(full(A+Dist(j)*U*S*V'))));
        Ranks(j)=size(S,1);
    end
    
    %% Gillis and Sharma method
    disp('Gillis-Sharma------------------------------')
    A=full(A);
    feval=@(x) 0.5*(sum(max(real(x+delta),0).^2));

    timmax=10;
    maxiter=1e16;
    algo=0;
    for j=totalsize-3:totalsize-1
        algo=algo+1;
        tic;
        [X,~,~]=neareststablealgo(A,maxiter,timmax,algo);
        CPUtimes(j)=toc;
        Delta_GS=X-A;
        Dist(j)=norm(Delta_GS,'fro');
        eig_tot_GS=eig(X);
        Fun(j)=feval(eig_tot_GS);
        Abscissa(j)=max(real(eig_tot_GS));
        Ranks(j)=rank(Delta_GS);
    end
    
    %% Noferini and Poloni method
    disp('Noferini-Poloni------------------------------')
    j=totalsize;
    tic;
    [B,~,~,~,T]=nearest_stable_complex(A,'hurwitz');
    CPUtimes(j)=toc;
    Delta_NP=B-A;
    Dist(j)=norm(Delta_NP,'fro');
    eig_tot_NP=diag(T);
    Fun(j)=feval(eig_tot_NP);
    Abscissa(j)=max(real(eig_tot_NP));
    Ranks(j)=rank(Delta_NP);
    
    %% Tables construction
    close all
    figure(1)
    
    rownames={[nex+4,1]};
    rownames{1}='Std Adaptive Rank';    
    for j=2:nex
        rk=r+j-2;
        rownames{j}=['Std Fixed Rank (' num2str(rk) ')'];
    end
    rownames{totalsize-3}='GS-BCD';
    rownames{totalsize-2}='GS-Grad';
    rownames{totalsize-1}='GS-FGM';
    rownames{totalsize}='NP';
    columnnames={'Distance','Rank','Functional','Abscissa','CPUtime'};
    
    Res={[5,1]};
    Res{1}=Dist;
    Res{2}=Ranks;
    Res{3}=Fun;
    Res{4}=Abscissa;
    Res{5}=CPUtimes;
    T1=table(Res{1},Res{2},Res{3},Res{4},Res{5},'RowNames',rownames); 
    UT1=uitable('Data',T1{:,:},'ColumnName',columnnames,...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

