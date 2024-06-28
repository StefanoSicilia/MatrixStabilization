%% Test for the Illustrative matrix
    
    %% Example 10-by-10 with 6 unstable eigenvalues
    n=10;
    r=6;
    A=[0,1,1,1,-1,0,-1,0,0,0;
        1,-1,0,1,1,0,1,0,0,0;
        -1,0,-1,-1,-1,1,1,1,0,0;
        1,0,0,-1,1,-1,-1,1,0,0;
        0,0,-1,1,0,1,1,-1,0,0;
        0,-1,1,1,-1,0,0,1,1,0;
        -1,1,-1,1,1,0,-1,0,1,1;
        0,0,1,-1,-1,1,1,1,-1,1;
        0,0,0,0,0,0,0,-1,1,-1;
        0,0,0,0,0,0,0,0,-1,1;];
    
    %% Parameters of the method
    epsilon=2;
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
    
    epsilon=2.4;
    for j=1:6
        parameters.ranktol=10^(-j);
        [~,~,~,F,~,Ranks]=InnerIter_UA(A,epsilon,parameters);
        Results(j,3:4)=[F(end),max(Ranks)];
    end 
    disp(Results)
    
    