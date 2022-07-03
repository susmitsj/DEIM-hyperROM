function mu = getmu(neq,r,ak,rk,Vr,Q,...
    uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0)
    
% Function to obtain the function in modified newton's method, mu(x) =
% f(x)/f'(x)
% Inputs: No. of equations neq
%         Reduced order dimension r
%           As required by the fdjac function


% Obtain the Jacobian J = f'(x)

%Determine the Jacobian of rk w.r.t. ak
J = fdjac(ak,rk,Vr,Q,...
    uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0);
    
mu = zeros(neq*r,1); % Initiate mu(x)
    
% Initiate equation index
k = 1;
while k <= 3 && norm(J((k-1)*r+1:k*r,(k-1)*r+1:k*r)-zeros(r,r),2)>1e-12
    if rank(J((k-1)*r+1:k*r,(k-1)*r+1:k*r)) ~= length(ak)
        mu((k-1)*r+1:k*r) = ...
        pinv(J((k-1)*r+1:k*r,(k-1)*r+1:k*r))*rk((k-1)*r+1:k*r);
    else
        mu((k-1)*r+1:k*r) = ...
        J((k-1)*r+1:k*r,(k-1)*r+1:k*r)\rk((k-1)*r+1:k*r);
    end
    % Update k
    k = k+1;
end

end

