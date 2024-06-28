%% Test for Grcar matrix
    
    %% Grcar 20-by-20 with 20 unstable eigenvalues
    n=20;
    r=n;
    A=gallery('grcar',n);
    v=sort(eig(A),'descend','ComparisonMethod','real');
    
    %% Parameters of the method
    epsilon=4;
    delta=1e-3;
    maxit=100;
    h=1;
    tol=1e-9;
    theta=1.3;
    safestop=15;
    parameters=struct('delta',delta,'maxit',maxit,'h',h,'tol',tol,...
        'theta',theta,'r',r,'safestop',safestop,'fun','standard');
    
    %% Main Computation
    Results=zeros(6,4);
    close all
    for j=1:6
        parameters.ranktol=10^(-j);
        [~,~,~,F,~,Ranks]=InnerIter_UA(A,epsilon,parameters);
        Results(j,1:2)=[F(end),max(Ranks)];
    end    
    
    epsilon=5.5;
    for j=1:6
        parameters.ranktol=10^(-j);
        [~,~,~,F,~,Ranks]=InnerIter_UA(A,epsilon,parameters);
        Results(j,3:4)=[F(end),max(Ranks)];
    end 
    disp(Results)  
    
    