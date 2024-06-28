function eta=ComputEta(L,Gamma,R,U,V,E)
%% ComputeEta:
% Computes the real part of the Frobenius inner product between P_Y G and 
% E, where Y=U*S*V' and G=L*Gamma*R'.
% Its the value eta for the low-rank-adaptove ODE.

    RGamma=R*Gamma';
    UE=U'*E;
    LU=L'*U;
    VRGamma=V'*RGamma;
    LEV=L'*E*V;
    eta=real(trace(UE*RGamma*LU+VRGamma*LEV-VRGamma*LU*UE*V));

end