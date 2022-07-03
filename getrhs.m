function rk = getrhs(uref,Vr,ak,Q,uinf, rho, imax, jmax, imms, isgs,...
    loc_p, loc_u, loc_v)

% Function to obtain the RHS of the governing equations starting from a
% reduced vector
% Inputs: uref(neq*Nx1) Reference solution vector
%         Vr (neq*N x neq*r) Block diagonal test basis matrix
%         ak (neq*rx1) Reduced solution vector
%         Q (neq*r x neq*N) Block diagonal interpolating matrix
%         Rest of the inputs as required by cavity_solver_iter
% Outputs: rk (neq*rx1) RHS vector

% No. of equations
neq = 3; %p, u, v
% Full order dimension
N = length(uref)/neq;

%Create u (neq*Nx1) from generalized co-ordinates
u = uref + Vr*ak;
% Set negative pressure values to zero
u(1:N)=max(1e-15,u(1:N));
% Separate the variables
u_mat = zeros(N, neq);

for s = 1 : neq
    u_mat(:,s) = u((s-1)*N+1 : s*N);
end
u_p = u_mat(:,1);
u_u = u_mat(:,2);
u_v = u_mat(:,3);

% Create a (imax x jmax x 3) matrix u
uin = zeros(imax,jmax,neq);

for s = 1 : jmax
    uin(:,s,1) = u_p((s-1)*imax+1:s*imax);
    uin(:,s,2) = u_u((s-1)*imax+1:s*imax);
    uin(:,s,3) = u_v((s-1)*imax+1:s*imax);
end
        
%Get RHS (neq*rx1)
[~, ~, ~, rhs_p, rhs_u, rhs_v, ~] =...
    cavity_solver_iter(uinf, rho, imax, jmax, imms, isgs, uin,...
    loc_p, loc_u, loc_v);
rk = Q*[rhs_p; rhs_u; rhs_v];

end