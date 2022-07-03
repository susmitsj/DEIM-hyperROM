function u0 = getic(x, u_s, Z, X)

%Function to obtain the initial condition (mean centered) from the snapshots
%via interpolation
%Inputs: New parameter vector x (uinf, mu)
%       Solution snapshot matrix, ps
%       Snapshot right singular vectors Z
%       Snapshot parameter values
%Outputs: Interpolated initial condition, u

tol = 1e-8;

% No. of snapshots
N_s = size(u_s,2); 

%Initiate variables
d=zeros(N_s,1);

for i = 1 : N_s
    d(i)=norm(x-X(i),2)^2;
end

if min(d)<tol
   u0 = u_s(:,i);
else
    a = 1./d;
    c = sum(a);
    a = a/c;
    u0 = u_s*a;
end

%Mean centered initial condition
uref = mean(u_s,2);

%Projected initial condition
u0 = uref + Z*Z'*(u0-uref);

