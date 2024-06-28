function [y,z]=Psi_HI(a,b,x)
%% Psi_HI: HermiteInterpolant
% Let p be the cubic polynomial such that p(a)=1, p(b)=p'(a)=p'(b)=0.
% Define          1   if x<a
%       psi(x) = p(x) if a<=x<=b
%                 0   if x>b
% Then y=psi(x) and z=psi'(x).
% It also accepts x as a vector, by acting entrywise.

    den=(a-b)^3;
    pindex=(x>a & x<b);
    xp=x(pindex);
    l=length(x);
    y=zeros(l,1);
    y(x>=b)=1;
    y(pindex)=(xp-a).^2.*(2*xp-3*b+a)/den;
    z=zeros(l,1);
    z(pindex)=6*((xp-b).*(xp-a))./den;

end