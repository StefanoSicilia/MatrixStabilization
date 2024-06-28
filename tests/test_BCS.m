%% Test for BCS matrix

    %% Matrix structure and eigenvalues
    n=13992;
    nnzhalved=316740;
    nnz=633480;
    r=2;
    A=load('bcsstk29.mtx');
    A=sparse(A(2:end,1),A(2:end,2),ones(nnzhalved,1),n,n);
    A=A+A'-diag(diag(A))-51.5*speye(n);
    v=sort(eigs(A,r+5,'largestreal'),'descend','ComparisonMethod','real');
    [row,col]=find(A);
    
    %% Parameters of the inner iteration
    delta=1e-3;
    maxit=150; 
    h=1e-4;
    tol=1e-9;
    theta=1.3;
    safestop=15;
    ranktol=1e-9;
    
    %% Parameters of the outer iteration
    el=1e-3;
    eu=5;
    e0=1e-1;
    epsilon=[el,eu,e0];
    niter=200; 
    toler=tol;   
    
    %% Grouping of the parameters
    InnerPar=struct('delta',delta,'maxit',maxit,'h',h,'tol',tol,...
        'theta',theta,'r',r,'safestop',safestop,'fun','standard',...
        'ranktol',ranktol,'row',row,'col',col);
    OuterPar=struct('epsilon',epsilon,'niter',niter,'toler',toler);
    
    %% Main Computation    
    % Adaptive
    OuterPar.type='SA';
    tic;
    [Dist,U,S,V,Fun]=Stabilize(A,InnerPar,OuterPar);
    CPUtimes=toc;
    E=projsparse(row,col,U,S,V);
    eig_tot=eigs(A+Dist*E,r+1,'largestreal');
    eig_stable=eig_tot(eig_tot<-delta);
    eig_unclear=eig_tot(eig_tot<=0 & eig_tot>=-delta);
    eig_unstable=eig_tot(eig_tot>0);
    Abscissa=max(real(eig_tot));
    Ranks=size(S,1);
    
    %% Tables construction
    close all
    figure(1)
    
    rownames={'Std Adaptive Rank'};    
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
    p{pcount}=plot(real(eig_stable),imag(eig_stable),...
        marker{j},'Color',[0,0.5,0],'LineWidth',lw);
    hold on
    p{pcount+1}=plot(real(eig_unclear),imag(eig_unclear),...
        marker{j},'Color',[1,0.75,0],'LineWidth',lw);
    hold on
    p{pcount+2}=plot(real(eig_unstable),imag(eig_unstable),...
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
    legend(q,legendlabel_eig(current),'interpreter','latex','FontSize',12)    
    
    