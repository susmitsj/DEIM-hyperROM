function [s,M]=q_deim(U)

% Input  : U (N x m) with orthonormal columns
% Output : S selection of m row indices with guaranteed upper bound
%           norm(inv(U(S,:))) <= sqrt(N-m+1) * O(2^m).
%        : M is the matrix U*inv(U(S,:))  (N x m);
% The Q-DEIM projection of an N x 1 vector f is M*f(S).
% Coded by Zlatko Drmac, April 2015.
[n,m] = size(U) ;
if nargout == 1
    [~,~,P] = qr(U','vector') ; s = P(1:m) ;
else
    [~,R,P] = qr(U','vector') ; s = P(1:m) ;
    M = [eye(m) ; (R(:,1:m)\R(:,m+1:n))'] ;
    Pinverse(P) = 1 : n ; M = M(Pinverse,:) ;
end
% I=eye(n);
% S=zeros(n,m);
% for j=1:m
%     S(:,j)=I(:,s(j));
% end
s=s';
end
