function [S,M,s]=deim(u)

% Input  : U (N x m) with orthonormal columns
% Output : S selection of m row indices with guaranteed upper bound
%           norm(inv(U(S,:))) <= sqrt(N-m+1) * O(2^m).
%        : M is the matrix U*inv(U(S,:))  (N x m);
% The Q-DEIM projection of an N x 1 vector f is M*f(S).
% Coded by Chaturantabut and Sorensen.

[N,m] = size(u) ;

%Initiate Variables
rho=zeros(m,1);
s=zeros(m,1);

[rho(1),s(1)]=max(abs(u(:,1)));
U=u(:,1);
I=eye(N);
S=I(:,s(1));

for j=2:m
    c=(S'*U)\(S'*u(:,j));
    r=u(:,j)-U*c;
    [rho(j),s(j)]=max(abs(r));
    U=[U,u(:,j)];
    S=[S,I(:,s(j))];
end

M=U*inv(S'*U);
