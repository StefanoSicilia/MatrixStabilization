%% Test for Fidap matrix

    %% Matrix structure and eigenvalues
    n=1601;
    r=4;
    A=load('fidap004.mtx');
    A=sparse(A(2:end,1),A(2:end,2),A(2:end,3),n,n);
    A=A-1.5*speye(n);
    v=sort(eig(full(A)),'descend','ComparisonMethod','real');
    [row,col]=find(A);
    
    %% Parameters of the inner iteration
    delta=1e-3;
    maxit=150; 
    h=1e-2;
    tol=1e-9;
    theta=1.3;
    safestop=15;
    stdfun='standard';
    Hermitefun='Hermite';
    type='SA';
    ranktol=1e-9;
    
    %% Parameters of the outer iteration
    el=1e-3;
    eu=0.5;
    e0=1e-1;
    epsilon=[el,eu,e0];
    niter=200; 
    toler=tol;   
    
    %% Grouping of the parameters
    InnerParsStd=struct('delta',delta,'maxit',maxit,'h',h,'tol',tol,...
        'theta',theta,'r',r,'safestop',safestop,'fun','standard',...
        'ranktol',ranktol,'row',row,'col',col);
    InnerParsHer=struct('delta',delta,'maxit',maxit,'h',h,'tol',tol,...
        'theta',theta,'safestop',safestop,'r',r,'fun','Hermite',...
        'ranktol',ranktol,'row',row,'col',col);
    OuterPar=struct('epsilon',epsilon,'niter',niter,'toler',toler);
    
    %% Structures for the outputs and plots 
    % Sizes of samples
    nex=3;
    totalsize=2*nex;
    
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
    OuterPar.type='SA';
    j=1;
    tic;
    [Dist(j),U,S,V,Fun(j)]=Stabilize(A,InnerParsStd,OuterPar);
    CPUtimes(j)=toc;
    E=projsparse(row,col,U,S,V);
    eig_tot{j}=eig(full(A+Dist(j)*E));
    eig_stable{j}=eig_tot{j}(eig_tot{j}<-delta);
    eig_unclear{j}=eig_tot{j}(eig_tot{j}<=0 & eig_tot{j}>=-delta);
    eig_unstable{j}=eig_tot{j}(eig_tot{j}>0);
    Abscissa(j)=max(real(eig_tot{j}));
    Ranks(j)=size(S,1);
    
    % Fixed rank
    OuterPar.type='SF';
    for j=2:nex
        rk=r+j-2;
        InnerParsStd.r=rk;
        tic;
        [Dist(j),U,S,V,Fun(j)]=Stabilize(A,InnerParsStd,OuterPar);
        CPUtimes(j)=toc;
        E=projsparse(row,col,U,S,V);
        Abscissa(j)=max(real(eig(full(A+Dist(j)*E))));
        Ranks(j)=size(S,1);
    end
    
    %% Main Computation: Hermite functional
    disp('Hermite-------------------------------')
    
    % Adaptive
    OuterPar.type='SA';
    j=nex+1;
    tic;
    [Dist(j),U,S,V,Fun(j)]=Stabilize(A,InnerParsHer,OuterPar);
    CPUtimes(j)=toc;
    E=projsparse(row,col,U,S,V);
    k=2;
    eig_tot{k}=eig(full(A+Dist(j)*E));
    eig_stable{k}=eig_tot{k}(eig_tot{k}<-delta);
    eig_unclear{k}=eig_tot{k}(eig_tot{k}<=0 & eig_tot{k}>=-delta);
    eig_unstable{k}=eig_tot{k}(eig_tot{k}>0);
    Abscissa(j)=max(real(eig_tot{k}));
    Ranks(j)=size(S,1);
    
    % Fixed rank
    OuterPar.type='SF';
    for j=(nex+2):totalsize
        rk=r+j-2-nex;
        InnerParsHer.r=rk;
        tic;
        [Dist(j),U,S,V,Fun(j)]=Stabilize(A,InnerParsHer,OuterPar);
        CPUtimes(j)=toc;
        E=projsparse(row,col,U,S,V);
        Abscissa(j)=max(real(eig(full(A+Dist(j)*E))));
        Ranks(j)=size(S,1);
    end
    
    %% Tables construction
    close all
    figure(1)
    
    rownames={[totalsize,1]};
    rownames{1}='Std Adaptive Rank';    
    for j=2:nex
        rk=r+j-2;
        rownames{j}=['Std Fixed Rank (' num2str(rk) ')'];
    end
    rownames{nex+1}='Her Adaptive Rank';
    for j=(2+nex):totalsize
        rk=r+j-2-nex;
        rownames{j}=['Her Fixed Rank (' num2str(rk) ')'];
    end
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
    
    %% Plot    
    % Plots parameters
    marker={'s'};
    nfig=1;
    trinfig=nfig*3+1;
    p={[trinfig,1]};
    legendlabel_eig={[trinfig,1]};
    legendlabel_eig{1}='Original';
    legendlabel_eig{2}='Adaptive Rank stable';
    legendlabel_eig{3}='Adaptive Rank unclear';
    legendlabel_eig{4}='Adaptive Rank unstable';
    lw=1.3;
    
    figure(2)
    p{1}=plot(real(v),imag(v),'ko','LineWidth',lw);
    j=1;
    pcount=3*j-1;
    hold on
    p{pcount}=plot(real(eig_stable{j}),imag(eig_stable{j}),...
        marker{j},'Color',[0,0.5,0],'LineWidth',lw);
    hold on
    p{pcount+1}=plot(real(eig_unclear{j}),imag(eig_unclear{j}),...
        marker{j},'Color',[1,0.75,0],'LineWidth',lw);
    hold on
    p{pcount+2}=plot(real(eig_unstable{j}),imag(eig_unstable{j}),...
        marker{j},'Color',[1,0,0],'LineWidth',lw);
    ax=gca;
    plot([-delta,-delta],[ax.YLim(1),ax.YLim(2)],'b-')  
    q=zeros(trinfig,1);
    k=1;
    current=zeros(trinfig,1);
    for j=1:trinfig
       if ~isempty(p{j}) 
           current(k)=j;
           q(k)=p{j};
           k=k+1;
       end          
    end
    current=current(1:k-1);
    q=q(1:k-1);
    legend(q,legendlabel_eig(current),'interpreter','latex')
    title('Standard')
    
    
    figure(3)
    p{1}=plot(real(v),imag(v),'ko','LineWidth',lw);
    j=1;
    pcount=3*j-1;
    k=j+nfig;
    hold on
    p{pcount}=plot(real(eig_stable{k}),imag(eig_stable{k}),...
        marker{j},'Color',[0,0.5,0],'LineWidth',lw);
    hold on
    p{pcount+1}=plot(real(eig_unclear{k}),imag(eig_unclear{k}),...
        marker{j},'Color',[1,0.75,0],'LineWidth',lw);
    hold on
    p{pcount+2}=plot(real(eig_unstable{k}),imag(eig_unstable{k}),...
        marker{j},'Color',[1,0,0],'LineWidth',lw);
    ax=gca;
    plot([-delta,-delta],[ax.YLim(1),ax.YLim(2)],'b-')  
    q=zeros(trinfig,1);
    k=1;
    current=zeros(trinfig,1);
    for j=1:trinfig
       if ~isempty(p{j}) 
           current(k)=j;
           q(k)=p{j};
           k=k+1;
       end          
    end
    current=current(1:k-1);
    q=q(1:k-1);
    legend(q,legendlabel_eig(current),'interpreter','latex')
    title('Hermite')
    
       
    

    
    
    
    
    
    
    
    
    
    