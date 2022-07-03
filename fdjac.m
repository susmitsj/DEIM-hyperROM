function [J] = fdjac(ak,rk,Vr,Q,...
    uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0)

%Function to calculate the finite difference Jacobian of residual w.r.t.
%the generalized co-ordinate vector
%Inputs: Small increment, e
%       Generalized co-ordinates ak (neq*rx1)
%       Reference state, ubar (neq*Nx1)
%       Residuals, rk (neq*rx1)
%       POD Basis, Vr (neq*Nxneq*r)
%       RHS POD Basis, Q (neq*rxneq*N)
%       Mean of snapshots, uref (neq*N)
%Outputs: Jacobian J (neq*rxneq*r)

%No. of equations
neq = 3;
%Reduced order dimension
r = length(ak)/neq;

%Initiate the Jacobian of rk w.r.t. ak 
%Block diagonal matrix
% J = eye(neq*r); % Initialized as identity matrix
J = zeros(neq*r,neq*r); % Initialized as zeros matrix

%Small change in ak
e = 1e-2*ak;

% Initiate the equation index
k = 1;
%Determine the Jacobian of rk w.r.t. ak
while k <= 3 && norm(rk((k-1)*r+1:k*r),2) > 1e-12
    for j = 1 : r
        %Get the corresponding ak to the equation
        akk = ak((k-1)*r+1:k*r);
        %Change each element by e for every j
        akk(j) = akk(j)+e((k-1)*r+j);
        %Create a (neq*rx1) from the new vector
        a = ak;
        a((k-1)*r+1:k*r) = akk;
        
        % Get RHS
        res = getrhs(uref,Vr,a,Q,uinf, rho, imax, jmax, imms, isgs,...
            loc_p, loc_u, loc_v);
        
        % Create Jacobian
        J((k-1)*r+1:k*r,(k-1)*r+j) = ...
            (res((k-1)*r+1:k*r)/d0(k)-rk((k-1)*r+1:k*r))/e((k-1)*r+j);
        
    end
    % Update equation index
    k = k+1;
end
    
end


