%% Test for the Smoke matrix
    
    %% Smoke n-by-n matrix
    n=20;
    r=n;
    A=gallery('smoke',n);
    v=sort(eig(A),'descend','ComparisonMethod','real');
    
    %% Parameters of the method
    epsilon=2.5;
    delta=1e-3;
    maxit=100;
    h=1;
    tol=1e-9;
    theta=1.3;
    safestop=15;
    parameters=struct('delta',delta,'maxit',maxit,'h',h,'tol',tol,...
        'theta',theta,'r',r,'safestop',safestop,'fun','standard');
    
    %% Main Computation
    Results=zeros(6,2);
    close all
    for j=1:6
        parameters.ranktol=10^(-j);
        [U,S,V,F,Times,Ranks]=InnerIter_UA(A,epsilon,parameters);
        if mod(j,2)==0
            figure
            semilogy(Times(2:end),F(2:end),'r-o','LineWidth',1)
            title('Objective functional versus time steps')
            figure
            plot(1:length(Ranks),Ranks,'b-s','LineWidth',1)
            curtick = get(gca, 'yTick');
            yticks(unique(round(curtick)));
            title('Ranks')
        end
        Results(j,:)=[F(end),max(Ranks)];
    end  
    disp(Results)
    
    