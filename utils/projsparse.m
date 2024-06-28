function A=projsparse(row,col,U,S,V)
%% projsparse:
% Projects the n-by-n rank-r matrix Y=U*S*V' onto the sparsity pattern 
% determined by row and col. The output A is a sparse real matrix.
% It assumes that U and V are n-by-r orthogonal matrices and that S is an
% r-by-r invertible matrix, with r<=n.


    %% Main computation
    n=size(U,1);
    VS=V*S;
    a=sum(((U(row,:).*conj(VS(col,:)))),2);
    A=real(sparse(row,col,a,n,n));
        
end