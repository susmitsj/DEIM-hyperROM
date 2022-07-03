function [rom_p, rom_u, rom_v, iter] = hrom(method, x, imax, jmax, imms, isgs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        DEIM-HROM Framework        %%%%%%%%%%%%
%%%%           Prepared by Susmit Joshi              %%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function to obtain solution by using Hyper-ROM
%Inputs: Method ('qdeim', 'deim')
%Outputs: Converged solution vector, u (neq*Nx1)
%       Iterations needed to solve the steady-state problem, iter

%% Load POD Basis

load('redbasis.mat','Vr_p','Vr_u','Vr_v',...
    'Ur_p','Ur_u','Ur_v');
load('snap.mat','X','ps','us','vs');
load('svd.mat','Z_p','Z_u','Z_v');

%Dimensions
N = size(Vr_p,1); %FOM dimension
r = size(Vr_u,2); %ROM dimension
m = size(Ur_v,2); %RHS interpolation dimension
neq = 3; 

%% Reference condition

%Find the mean
ps_m = mean(ps,2);
us_m = mean(us,2);
vs_m = mean(vs,2);

uref = [ps_m; us_m; vs_m];

%% Create a global Vr matrix

Vr = zeros(neq*N,neq*r);
Vr(1:N,1:r) = Vr_p;
Vr(N+1:2*N,r+1:2*r) = Vr_u;
Vr(2*N+1:3*N,2*r+1:3*r) = Vr_v;

%% Apply DEIM on RHS POD Basis

switch method
    case 'qdeim'
        [p_p,M_p] = q_deim(Ur_p);
        [p_u,M_u] = q_deim(Ur_u);
        [p_v,M_v] = q_deim(Ur_v);
        
    case 'deim'
        [~,M_p,p_p] = deim(Ur_p);
        [~,M_u,p_u] = deim(Ur_u);
        [~,M_v,p_v] = deim(Ur_v);
        
    case 'no'
        M_p = eye(N);
        M_u = eye(N);
        M_v = eye(N);
        
end

%% Precomputing RHS Multiplier

%Create a global Q matrix
Q = zeros(neq*r,neq*m);
Q(1:r,1:m) = Vr_p'*M_p;
Q(r+1:2*r,m+1:2*m) = Vr_u'*M_u;
Q(2*r+1:3*r,2*m+1:3*m) = Vr_v'*M_v;

%% Get initial conditions

% Get the initial condition
u_p0 = getic(x, ps, Z_p, X);  
u_u0 = getic(x, us, Z_u, X);    
u_v0 = getic(x, vs, Z_v, X);    

% Create a global u
u0 = [u_p0; u_u0; u_v0];

%Create a global a (neq*rx1)
a0 = Vr'*(u0-[mean(ps,2);mean(us,2);mean(vs,2)]);

% Obtain the RHS interpolation locations 

[loc_p] = rhsinter(jmax, p_p);
[loc_u] = rhsinter(jmax, p_u);
[loc_v] = rhsinter(jmax, p_v);

%Get RHS (neq*rx1)
uinf = x(1);
rho = x(2);

% Get RHS
r0 = getrhs(uref,Vr,a0,Q,uinf, rho, imax, jmax, imms, isgs,...
    loc_p, loc_u, loc_v);

%% Find the generalized co-ordinates such that projected RHS is 0

%Set convergence tolerance
tol = 1e-4;
%Set small difference for calculating the finite difference Jacobian (1% of
%the values in the generalized co-ordinates vector)
% e = 1e-4*a0;
%Secant method
[a, iter] = secant(tol,a0,r0,Vr,Q,0,...
    uinf, rho, imax, jmax, imms, isgs, uref,...
    loc_p, loc_u, loc_v);
% iter=1;
%% Find the converged solutions
u = uref + Vr*a;

% Separate the variables
u_mat = zeros(N, neq);

for j = 1 : neq
    u_mat(:,j) = u((j-1)*N+1 : j*N);
end
rom_p = u_mat(:,1);
rom_u = u_mat(:,2);
rom_v = u_mat(:,3);





