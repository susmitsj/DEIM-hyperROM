function [sol_p, sol_u, sol_v, rhs_p, rhs_u, rhs_v, p, U, V] = cavity_solver_mms(uinf, rho, imax, jmax, isgs, sim_no)
tic   %begin timer function
%--- Variables for file handling ---
%--- All files are globally accessible ---

global fp1 % For output of iterative residual history
global fp2 % For output of field data (solution)
%   global fp3 % For writing the restart file
%   global fp4 % For reading the restart file
%   global fp5 % For output of final DE norms (only for MMS)
%$$$$$$   global fp6 % For debug: Uncomment for debugging.

global neq nmax imms
global zero tenth sixth fifth fourth third half one two three four six
global iterout irstr ipgorder lim cfl Cx Cy toler rkappa Re pinf rhoinv xmin xmax ymin ymax Cx2 Cy2 fsmall
global rlength rmu vel2ref dx dy rpi phi0 phix phiy phixy apx apy apxy fsinx fsiny fsinxy
%--------------------------------------------------------------------------------------------------------------------------%
global rhs
%--------------------------------------------------------------------------------------------------------------------------%

%**Use these variables cautiously as these are globally accessible from all functions.**

global u;         % Solution vector [p, u, v]^T at each node
global uold;      % Previous (old) solution vector
global s;         % Source term
global dt;        % Local time step at each node
global artviscx;  % Artificial viscosity in x-direction
global artviscy;  % Artificial viscosity in y-direction
global ummsArray; % Array of umms values (funtion umms evaluated at all nodes)

%************ Following are fixed parameters for array sizes *************
%--------------------------------------------------------------------------------------------------------------------------%
% imax = 65;   % Number of points in the x-direction (use odd numbers only)
% jmax = 65;   % Number of points in the y-direction (use odd numbers only)
%--------------------------------------------------------------------------------------------------------------------------%
neq = 3;     % Number of equation to be solved ( = 3: mass, x-mtm, y-mtm)
%********************************************
%***** All  variables declared here. **
%**** These variables SHOULD not be changed *
%********* by the program once set. *********
%********************************************
%**** The variables declared "" CAN ****
%** not be changed by the program once set **
%********************************************

%--------- Numerical constants --------
zero   = 0.0;
tenth  = 0.1;
sixth  = 1.0/6.0;
fifth  = 0.2;
fourth = 0.25;
third  = 1.0/3.0;
half   = 0.5;
one    = 1.0;
two    = 2.0;
three  = 3.0;
four   = 4.0;
six    = 6.0;

%--------- User sets inputs here  --------

nmax = 500000;        % Maximum number of iterations
iterout = 5000;       % Number of time steps between solution output
%--------------------------------------------------------------------------------------------------------------------------%
imms = 1;             % Manufactured solution flag: = 1 for manuf. sol., = 0 otherwise
% isgs = 1;             % Symmetric Gauss-Seidel  flag: = 1 for SGS, = 0 for point Jacobi
%--------------------------------------------------------------------------------------------------------------------------%
irstr = 0;            % Restart flag: = 1 for restart (file 'restart.in', = 0 for initial run
ipgorder = 0;         % Order of pressure gradient: 0 = 2nd, 1 = 3rd (not needed)
lim = 1;              % variable to be used as the limiter sensor (= 1 for pressure)

cfl  = 0.5;      % CFL number used to determine time step
Cx = 0.01;     	% Parameter for 4th order artificial viscosity in x
Cy = 0.01;      	% Parameter for 4th order artificial viscosity in y
toler = 1.e-4; 	% Tolerance for iterative residual convergence
rkappa = 0.1;   	% Time derivative preconditioning constant
pinf = 0.801333844662; % Initial pressure (N/m^2) -> from MMS value at cavity center
%--------------------------------------------------------------------------------------------------------------------------%
% uinf = 1.0;               % Lid velocity (m/s)
% rho = 1.0;                % Density (kg/m^3)
rmu = 5e-4;              % Constant viscosity (N*s/m^2)
L = 0.05;                   % Maximum length (m)
Re = rho*uinf*L/rmu;        % Reynolds number = rho*Uinf*L/rmu
ymax = L;                   % Maximum y location (m)
xmax = L;                   % Maximum x location (m)
%--------------------------------------------------------------------------------------------------------------------------%
xmin = 0.0;      % Cavity dimensions...: minimum x location (m)
ymin = 0.0;      %                       maximum y location (m)
Cx2 = 0.0;       % Coefficient for 2nd order damping (not required)
Cy2 = 0.0;     	 % Coefficient for 2nd order damping (not required)
fsmall = 1.e-20; % small parameter

%-- Derived input quantities (set by function 'set_derived_inputs' called from main)----

rhoinv =  -99.9; 	% Inverse density, 1/rho (m^3/kg)
rlength = -99.9;  	% Characteristic length (m) [cavity width]
vel2ref = -99.9;  	% Reference velocity squared (m^2/s^2)
dx = -99.9; 		% Delta x (m)
dy = -99.9;  		% Delta y (m)
rpi = -99.9; 		% Pi = 3.14159... (defined below)

%-- constants for manufactured solutions ----
phi0 = [0.25, 0.3, 0.2];          % MMS constant
phix = [0.5, 0.15, 1.0/6.0];      % MMS amplitude constant
phiy = [0.4, 0.2, 0.25];          % MMS amplitude constant
phixy = [1.0/3.0, 0.25, 0.1];     % MMS amplitude constant
apx = [0.5, 1.0/3.0, 7.0/17.0]; 	% MMS frequency constant
apy = [0.2, 0.25, 1.0/6.0];         % MMS frequency constant
apxy = [2.0/7.0, 0.4, 1.0/3.0];     % MMS frequency constant
fsinx = [0.0, 1.0, 0.0];            % MMS constant to determine sine vs. cosine
fsiny = [1.0, 0.0, 0.0];            % MMS constant to determine sine vs. cosine
fsinxy = [1.0, 1.0, 0.0];           % MMS constant to determine sine vs. cosine
% Note: fsin = 1 means the sine function
% Note: fsin = 0 means the cosine function
% Note: arrays here refer to the 3 variables

%************************************************************************
%      						Main Function
%************************************************************************
%----- Looping indices --------

i = 0;                       % i index (x direction)
j = 0;                       % j index (y direction)
k = 0;                       % k index (# of equations)
n = 0;	                   % Iteration number index

conv = -99.9 ; % Minimum of iterative residual norms from three equations


%--------- Solution variables declaration --------

ninit = 0;        	% Initial iteration number (used for restart file)

%$$$$$$    u(imax,jmax,neq);         % Solution vector (p, u, v)^T at each node
%$$$$$$    uold(imax,jmax,neq);      % Previous (old) solution vector
%$$$$$$    s(imax,jmax,neq);    % Source term
%$$$$$$    dt(imax,jmax);        % Local time step at each node
%$$$$$$    artviscx(imax,jmax);  % Artificial viscosity in x-direction
%$$$$$$    artviscy(imax,jmax);  % Artificial viscosity in y-direction
res = [0,0,0];              % Iterative residual for each equation
resinit = [0,0,0];          % Initial iterative residual for each equation (from iteration 1)
rL1norm = [0,0,0];          % L1 norm of discretization error for each equation
rL2norm = [0,0,0];          % L2 norm of discretization error for each equation
rLinfnorm = [0,0,0];        % Linfinity norm of discretization error for each equation
rtime = -99.9;         % Variable to estimate simulation time
dtmin = 1.0e99;        % Minimum time step for a given iteration (initialized large)

x = -99.9;       % Temporary variable for x location
y = -99.9;       % Temporary variable for y location

% Solution variables initialization with dummy values
% for i=1:imax
%  for j=1:jmax
%      dt(i,j) = -99.9;
%      artviscx(i,j) = -99.9;
%      artviscy(i,j) = -99.9;
%    for k=1:neq
%      u(i,j,k) = -99.9;
%      uold(i,j,k) = -99.9;
%      s(i,j,k) = -99.9;
%      res(k) = -99.9;
%      resinit(k) = -99.9;
%      res(k) = -99.9;
%      rL1norm(k) = -99.9;
%      rL2norm(k) = -99.9;
%      rLinfnorm(k) = -99.9;
%    end
%  end
% end

dt = zeros(imax,jmax);
artviscx = zeros(imax,jmax);
artviscy = zeros(imax,jmax);
u = zeros(imax,jmax,neq);
%--------------------------------------------------------------------------------------------------------------------------%
rhs = zeros(imax,jmax,neq);
%--------------------------------------------------------------------------------------------------------------------------%
uold = zeros(imax,jmax,neq);
ummsArray = zeros(imax,jmax,neq);
s = zeros(imax,jmax,neq);
res = zeros(neq,1);
resinit = zeros(neq,1);
rL1norm = zeros(neq,1);
rL2norm = zeros(neq,1);
rLinfnorm = zeros(neq,1);

dt(:,:) = -99.9;
artviscx(:,:) = -99.9;
artviscy(:,:) = -99.9;
u(:,:,:) = -99.9;
uold(:,:,:) = -99.9;
s(:,:,:) = -99.9;
res(:) = -99.9;
resinit(:) = -99.9;
rL1norm(:) = -99.9;
rL2norm(:) = -99.9;
rLinfnorm(:) = -99.9;


% Debug output: Uncomment and modify if debugging
%$$$$$$ fp6 = fopen("./Debug.dat","w");
%$$$$$$ fprintf(fp6,"TITLE = \"Debug Data Data\"\n");
%$$$$$$ fprintf(fp6,"variables=\"x(m)\"\"y(m)\"\"visc-x\"\"visc-y\"\n");
%$$$$$$ fprintf(fp6, "zone T=\"n=%d\"\n",n);
%$$$$$$ fprintf(fp6, "I= %d J= %d\n",imax, jmax);
%$$$$$$ fprintf(fp6, "DATAPACKING=POINT\n");

% Set derived input quantities
set_derived_inputs(uinf, rho, imax, jmax);

% Set up headers for output files
output_file_headers(imms, sim_no);

% Set Initial Profile for u vector
[ninit, rtime, resinit] = initial(ninit, rtime, resinit, imax, jmax, uinf);

% Set Boundary Conditions for u
set_boundary_conditions(imms, imax, jmax, uinf)

% Write out inital conditions to solution file
write_output( n,  resinit,  rtime, imax, jmax, imms);

% Initialize Artificial Viscosity arrays to zero (note: artviscx(i,j) and artviscy(i,j)
artviscx(:,:) = zero;
artviscy(:,:) = zero;

% Evaluate Source Terms Once at Beginning
%(only interior points; will be zero for standard cavity)
compute_source_terms(imax, jmax, imms, rho);

%========== Main Loop ==========
isConverged = 0;

for n = ninit:nmax
    % Calculate time step
    dtmin = compute_time_step(dtmin, imax, jmax, rho);
    
    % Save u values at time level n (u and uold are 2D arrays)
    uold = u;
    
    if isgs==1 % ==Symmetric Gauss Seidel==
        
        % Artificial Viscosity
        Compute_Artificial_Viscosity(imax, jmax, rho);
        
        % Symmetric Gauss-Siedel: Forward Sweep
        SGS_forward_sweep(imax, jmax, rho);
        
        % Set Boundary Conditions for u
        set_boundary_conditions(imms, imax, jmax, uinf);
        
        % Artificial Viscosity
        Compute_Artificial_Viscosity(imax, jmax, rho);
        
        % Symmetric Gauss-Siedel: Backward Sweep
        SGS_backward_sweep(imax, jmax, rho);
        
        % Set Boundary Conditions for u
        set_boundary_conditions(imms, imax, jmax, uinf);
    else
        if isgs==0 % ==Point Jacobi==
            
            % Artificial Viscosity
            Compute_Artificial_Viscosity(imax, jmax, rho);
            
            % Point Jacobi: Forward Sweep
            point_Jacobi(imax, jmax, rho);
            
            % Set Boundary Conditions for u
            set_boundary_conditions(imms, imax, jmax, uinf);
        else
            fprintf('ERROR: isgs must equal 0 or 1!\n');
            return;
        end
    end
    
    % Pressure Rescaling (based on center point)
    pressure_rescaling(imax, jmax, imms);
    
    % Update the time
    rtime = rtime + dtmin;
    
    % Check iterative convergence using L2 norms of iterative residuals
    [res, resinit, conv] = check_iterative_convergence(n, res, resinit, ninit, rtime, dtmin, imax, jmax);
    
    if(conv<toler)
        fprintf(fp1, '%d %e %e %e %e\n',n, rtime, res(1), res(2), res(3));
        isConverged = 1;
        break;
    end
    
    % Output solution and restart file every 'iterout' steps
    if( (mod(n,iterout)==0) )
        write_output( n,  resinit,  rtime, imax, jmax, imms);
    end
    
end  % ========== End Main Loop ==========

if isConverged == 0
    fprintf('Solution failed to converge in %d iterations!!!', nmax);
end

if isConverged == 1
    fprintf('Solution converged in %d iterations!!!', n);
end

% Calculate and Write Out Discretization Error Norms (will do this for MMS only)
Discretization_Error_Norms(rL1norm, rL2norm, rLinfnorm, imax, jmax, imms);

% Output solution and restart file
write_output( n,  resinit,  rtime, imax, jmax, imms);

% Close open files
fclose(fp1);
fclose(fp2);
%$$$$$$   fclose(fp6); % Uncomment for debug output (

%---------------------------------------------------------------------------------------------------------------------------%
% Output vectors
sol_p = zeros(imax*jmax,1);
sol_u = zeros(imax*jmax,1);
sol_v = zeros(imax*jmax,1);
rhs_p = zeros(imax*jmax,1);
rhs_u = zeros(imax*jmax,1);
rhs_v = zeros(imax*jmax,1);
p = zeros(imax*jmax,1);
U = zeros(imax*jmax,1);
V = zeros(imax*jmax,1);

for j = 1 : jmax
    sol_p((j-1)*imax + 1 : j*imax) = u(:,j,1);
    sol_u((j-1)*imax + 1 : j*imax) = u(:,j,2);
    sol_v((j-1)*imax + 1 : j*imax) = u(:,j,3);
    rhs_p((j-1)*imax + 1 : j*imax) = rhs(:,j,1);
    rhs_u((j-1)*imax + 1 : j*imax) = rhs(:,j,2);
    rhs_v((j-1)*imax + 1 : j*imax) = rhs(:,j,3);
    p((j-1)*imax + 1 : j*imax) = ummsArray(:,j,1);
    U((j-1)*imax + 1 : j*imax) = ummsArray(:,j,2);
    V((j-1)*imax + 1 : j*imax) = ummsArray(:,j,3);
end

%--------------------------------------------------------------------------------------------------------------------------%

% PrsMatrix = u(:,:,1);    %output arrays
% uvelMatrix = u(:,:,2);
% vvelMatrix = u(:,:,3);

toc  %end timer function
end

%**************************************************************************/
%*      					All Other	Functions					      */
%**************************************************************************/

%**************************************************************************
%**************************************************************************
function set_derived_inputs(uinf, rho, imax, jmax)
global one
global Re rhoinv xmin xmax ymin ymax
global rlength rmu vel2ref dx dy rpi

rhoinv = 1/rho;                            % Inverse density, 1/rho (m^3/kg) */
rlength = xmax - xmin;                       % Characteristic length (m) [cavity width] */
%--------------------------------------------------------------------------------------------------------------------------%
rmu = 5e-4;                               % Constant viscosity (N*s/m^2) 
% rmu = rho*uinf*rlength/Re;                   % Viscosity (N*s/m^2) */
%--------------------------------------------------------------------------------------------------------------------------%
vel2ref = uinf*uinf;                         % Reference velocity squared (m^2/s^2) */
dx = (xmax - xmin)/(imax - 1);          % Delta x (m) */
dy = (ymax - ymin)/(jmax - 1);          % Delta y (m) */
rpi = acos(-one);                            % Pi = 3.14159... */
fprintf('rho,V,L,mu,Re: %f %f %f %f %f\n',rho,uinf,rlength,rmu,Re);
end

%************************************************************************
function output_file_headers(imms, sim_no)

% Uses global variable(s): imms, fp1, fp2
% Note: The vector of primitive variables is:
%               u = [p, u, v]^T
% Set up output files (history and solution)

global fp1 fp2

%--------------------------------------------------------------------------------------------------------------------------%
% Set working directory
% d = 'C:\Users\susmi\Google Drive\Spring 2022\Academics\AOE 6444 V & V in SC\Project\Snapshots'; 

filename1 = sprintf('history_%d.dat,',sim_no);
fp1 = fopen(filename1,'w');
%--------------------------------------------------------------------------------------------------------------------------%
fprintf(fp1,'TITLE = "Cavity Iterative Residual History"\n');
fprintf(fp1,'variables="Iteration""Time(s)""Res1""Res2""Res3"\n');

%--------------------------------------------------------------------------------------------------------------------------%
filename2 = sprintf('cavity_%d.dat,',sim_no);
fp2 = fopen(filename2,'w');
%--------------------------------------------------------------------------------------------------------------------------%
fprintf(fp2,'TITLE = "Cavity Field Data"\n');
if (imms==1)
    
    fprintf(fp2,'variables="x(m)""y(m)""p(N/m^2)""u(m/s)""v(m/s)"');
    fprintf(fp2,'"p-exact""u-exact""v-exact""DE-p""DE-u""DE-v"\n');
    
else
    
    if (imms==0)
        
        fprintf(fp2,'variables="x(m)""y(m)""p(N/m^2)""u(m/s)""v(m/s)"\n');
        
    else
        
        fprintf('ERROR! imms must equal 0 or 1!!!\n');
        return;
    end
end

% Header for Screen Output
fprintf('Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum\n');

end
%************************************************************************
function [ninit, rtime, resinit] = initial(ninit, rtime, resinit, imax, jmax, uinf)
%
%Uses global variable(s): zero, one, irstr, imax, jmax, neq, uinf, pinf
%To modify: ninit, rtime, resinit, u, s

% i                        % i index (x direction)
% j                        % j index (y direction)
% k                        % k index (# of equations)
% x        % Temporary variable for x location
% y        % Temporary variable for y location

% This subroutine sets inital conditions in the cavity
% Note: The vector of primitive variables is:
%              u = (p, u, v)^T

global zero one irstr neq pinf xmax xmin ymax ymin
global u s ummsArray

% Initialize the ummsArray with values computed with umms function
for j=1:jmax
    for i=1:imax
        for k=1:neq
            x = (xmax - xmin)*(i-1)/(imax - 1);
            y = (ymax - ymin)*(j-1)/(jmax - 1);
            ummsArray(i,j,k) = umms(x,y,k);
        end
    end
end

if (irstr==0)   % Starting run from scratch
    ninit = 1;          % set initial iteration to one
    rtime = zero;       % set initial time to zero
    for k = 1:neq
        resinit(k) = one;
    end
    for j = 1:jmax
        for i = 1:imax
            u(i,j,1) = pinf;
            u(i,j,2) = zero;
            u(i,j,3) = zero;
            s(i,j,1) = zero;
            s(i,j,2) = zero;
            s(i,j,3) = zero;
        end
        u(i,jmax,2) = uinf; % Initialize lid (top) to freestream velocity
    end
else
    if (irstr==1)  % Restarting from previous run (file 'restart.in')
        fp4 = fopen('./restart.out','r'); % Note: 'restart.in' must exist!
        if (fp4==-1)
            fprintf('Error opening restart file. Stopping.\n');
            return;
        end
        fscanf(fp4, '%d %lf', ninit, rtime); % Need to known current iteration # and time value
        fscanf(fp4, '%lf %lf %lf', resinit(0), resinit(1), resinit(2)); % Needs initial iterative residuals for scaling
        for j=1:jmax
            for i=1:imax
                fscanf(fp4, '%lf %lf %lf %lf %lf', x, y, u(i,j,1), u(i,j,2), u(i,j,3));
            end
        end
        ninit = ninit + 1;
        fprintf('Restarting at iteration %d\n', ninit);
        fclose(fp4);
    else
        printf('ERROR: irstr must equal 0 or 1!\n');
        return;
    end
end

end
%************************************************************************
function set_boundary_conditions(imms, imax, jmax, uinf)
%
%Uses global variable(s): imms
%To modify: u (via other functions: bndry() and bndrymms())

% This subroutine determines the appropriate BC routines to call
if (imms==0)
    bndry(imax, jmax, uinf);
else
    if (imms==1)
        bndrymms(imax, jmax);
    else
        printf('ERROR: imms must equal 0 or 1!\n');
        return;
    end
end
end
%************************************************************************
function bndry(imax, jmax, uinf)
%
%Uses global variable(s): zero, one, two, half, imax, jmax, uinf
%To modify: u

% i                        % i index (x direction)
% j                        % j index (y direction)

global zero two half
global u

% This applies the cavity boundary conditions
% Side Walls
for j = 2:jmax-1
    %    u(1,j,1) = ( 18*u(2,j,1) - 9*u(3,j,1) + two*u(4,j,1) ) / 11;   % 3rd Order BC
    u(1,j,1) = two*u(2,j,1) - u(3,j,1);   % 2nd Order BC
    %    u(1,j,1) = u(2,j,1);      % 1st Order BC
    u(1,j,2) = zero;
    u(1,j,3) = zero;
    %   u(imax,j,1) = ( 18*u(imax-1,j,1) - 9*u(imax-2,j,1) \
    %                   + two*u(imax-3,j,1) ) / 11;  % 3rd Order BC
    u(imax,j,1) = two*u(imax-1,j,1) - u(imax-2,j,1);   % 2nd Order BC
    %  u(imax,j,1) = u(imax-1,j,1);    % 1st Order BC
    u(imax,j,2) = zero;
    u(imax,j,3) = zero;
    %$$$$$$     printf("j, u1, u2, u3: %d %f %f %f %f %f %f\n",j, u(0,j,0), u(0,j,1), u(0,j,2), \
    %$$$$$$     u(imax-1,j,0), u(imax-1,j,1), u(imax-1,j,2));
end

% Top/Bottom Walls
for i = 2:imax-1
    %   u(i,1,1) = ( 18*u(i,2,1) - 9*u(i,3,1) + two*u(i,4,1) ) / 11;   % 3rd Order BC
    u(i,1,1) = two*u(i,2,1) - u(i,3,1);   % 2nd Order BC
    %    u(i,1,1) = u(i,2,1);       % 1st Order BC
    u(i,1,2) = zero;
    u(i,1,3) = zero;
    %    u(i,jmax,1) = ( 18*u(i,jmax-1,1) - 9*u(i,jmax-2,1) \
    %                   + two*u(i,jmax-3,1) ) / 11;  % 3rd Order BC
    u(i,jmax,1) = two*u(i,jmax-1,1) - u(i,jmax-2,1);  % 2nd Order BC
    %    u(i,jmax,1) = u(i,jmax-1,1);     % 1st Order BC
    u(i,jmax,2) = uinf;
    u(i,jmax,3) = zero;
    %$$$$$$     printf("i, u1, u2, u3: %d %f %f %f %f %f %f\n",i, u(i,0,0), u(i,0,1), u(i,0,2), \
    %$$$$$$     u(i,jmax-1,0), u(i,jmax-1,1), u(i,jmax-1,2));
end

% Corners
u(1,jmax,1) = half*( u(2,jmax,1) + u(1,jmax-1,1) );
u(imax,jmax,1) = half*( u(imax-1,jmax,1) + u(imax,jmax-1,1) );
u(1,jmax,2) = zero;      % uinf   /two     uinf
u(imax,jmax,2) = zero;   % uinf   /two     uinf
u(1,jmax,3) = zero;
u(imax,jmax,3) = zero;

u(1,1,1) = half*( u(2,1,1) + u(1,2,1) );
u(imax,1,1) = half*( u(imax-1,1,1) + u(imax,2,1) );
u(1,1,2) = zero;
u(imax,1,2) = zero;
u(1,1,3) = zero;
u(imax,1,3) = zero;
end
%************************************************************************
function bndrymms(imax, jmax)
%
%Uses global variable(s): two, imax, jmax, neq, xmax, xmin, ymax, ymin, rlength
%To modify: u
% i                        % i index (x direction)
% j                        % j index (y direction)
% k                        % k index (# of equations)
% x        % Temporary variable for x location
% y        % Temporary variable for y location
% This applies the cavity boundary conditions for the manufactured solution

global two neq
global u ummsArray

% Side Walls
for j = 2:jmax-1
    i = 1;
    for k = 1:neq
        u(i,j,k) = ummsArray(i,j,k);
    end
    u(1,j,1) = two*u(2,j,1) - u(3,j,1);    % 2nd Order BC
    %    u(1,j,1) = u(2,j,1);                  % 1st Order BC
    
    i=imax;
    for k = 1:neq
        u(i,j,k) = ummsArray(i,j,k);
    end
    u(imax,j,1) = two*u(imax-1,j,1) - u(imax-2,j,1);   % 2nd Order BC
    %	u(imax,j,1) = u(imax-1,j,1);                       % 1st Order BC
end

% Top/Bottom Walls
for i=1:imax
    j = 1;
    for k = 1:neq
        u(i,j,k) = ummsArray(i,j,k);
    end
    u(i,1,1) = two*u(i,2,1) - u(i,3,1);   % 2nd Order BC
    %$$$$$$     u(i,1,1) = u(i,2,1);            % 1st Order BC
    
    j = jmax;
    for k = 1:neq
        u(i,j,k) = ummsArray(i,j,k);
    end
    u(i,jmax,1) = two*u(i,jmax-1,1) - u(i,jmax-2,1);   % 2nd Order BC
    %$$$$$$     u(i,jmax,1) = u(i,jmax-1,1);              % 1st Order BC
end
end
%************************************************************************
function [ummstmp] = umms(x, y, k)
%
%Uses global variable(s): one, rpi, rlength
%Inputs: x, y, k
%To modify: <none>
%Returns: umms

% ummstmp; % Define return value for umms as % precision

% termx       % Temp variable
% termy       % Temp variable
% termxy      % Temp variable
% argx        % Temp variable
% argy        % Temp variable
% argxy       % Temp variable

% This function returns the MMS exact solution

global one rpi rlength
global phi0 phix phiy phixy apx apy apxy fsinx fsiny fsinxy

argx = apx(k)*rpi*x/rlength;
argy = apy(k)*rpi*y/rlength;
argxy = apxy(k)*rpi*x*y/rlength/rlength;
termx = phix(k)*(fsinx(k)*sin(argx)+(one-fsinx(k))*cos(argx));
termy = phiy(k)*(fsiny(k)*sin(argy)+(one-fsiny(k))*cos(argy));
termxy = phixy(k)*(fsinxy(k)*sin(argxy)+(one-fsinxy(k))*cos(argxy));

ummstmp = phi0(k) + termx + termy + termxy;
end
%************************************************************************
function write_output(n,  resinit,  rtime, imax, jmax, imms)
%
%Uses global variable(s): imax, jmax, new, xmax, xmin, ymax, ymin, rlength, imms
%Uses global variable(s): ninit, u, dt, resinit, rtime
%To modify: <none>
%Writes output and restart files.

% i                        % i index (x direction)
% j                        % j index (y direction)
% k                        % k index (# of equations)

% x        % Temporary variable for x location
% y        % Temporary variable for y location

global xmax xmin ymax ymin
global u ummsArray
global fp2 fp3

% Field output
fprintf(fp2, 'zone T="n=%d"\n',n);
fprintf(fp2, 'I= %d J= %d\n',imax, jmax);
fprintf(fp2, 'DATAPACKING=POINT\n');

if (imms==1)
    for j=1:jmax
        for i=1:imax
            x = (xmax - xmin)*(i-1)/(imax - 1);
            y = (ymax - ymin)*(j-1)/(jmax - 1);
            fprintf(fp2,'%e %e %e %e %e %e %e %e %e %e %e\n', x, y, ...
                u(i,j,1), u(i,j,2), u(i,j,3), ummsArray(i,j,1), ummsArray(i,j,2), ummsArray(i,j,3), ...
                (u(i,j,1)-ummsArray(i,j,1)), (u(i,j,2)-ummsArray(i,j,2)), (u(i,j,3)-ummsArray(i,j,3)));
        end
    end
else
    if (imms==0)
        for j=1:jmax
            for i=1:imax
                x = (xmax - xmin)*(i-1)/(imax - 1);
                y = (ymax - ymin)*(j-1)/(jmax - 1);
                fprintf(fp2,'%e %e %e %e %e\n', x, y, ...
                    u(i,j,1), u(i,j,2), u(i,j,3));
            end
        end
    else
        fprintf('ERROR: imms must equal 0 or 1!\n');
        return;
    end
end

% Restart file: overwrites every 'iterout' iteration
fp3 = fopen('./restart.out','w');
fprintf(fp3,'%d %e\n', n, rtime);
fprintf(fp3,'%e %e %e\n', resinit(1), resinit(2), resinit(3));
for j=1:jmax
    for i=1:imax
        x = (xmax - xmin)*(i-1)/(imax - 1);
        y = (ymax - ymin)*(j-1)/(jmax - 1);
        fprintf(fp3,'%e %e %e %e %e\n', x, y, ...
            u(i,j,1), u(i,j,2), u(i,j,3));
    end
end
fclose(fp3);
end
%************************************************************************
function compute_source_terms(imax, jmax, imms, rho)
%
%Uses global variable(s): imax, jmax, imms, rlength, xmax, xmin, ymax, ymin
%To modify: s (source terms)

% i                        % i index (x direction)
% j                        % j index (y direction)

% x        % Temporary variable for x location
% y        % Temporary variable for y location

% Evaluate Source Terms Once at Beginning (only %erior po%s; will be zero for standard cavity)

global xmax xmin ymax ymin
global s

for j=2:jmax-1
    for i=2:imax-1
        x = (xmax - xmin)*(i-1)/(imax - 1);
        y = (ymax - ymin)*(j-1)/(jmax - 1);
        s(i,j,1) = (imms)*srcmms_mass(x, y, rho);
        s(i,j,2) = (imms)*srcmms_xmtm(x, y, rho);
        s(i,j,3) = (imms)*srcmms_ymtm(x, y, rho);
    end
end
end
%************************************************************************
function [srcmasstmp] = srcmms_mass(x, y, rho)
%
%Uses global variable(s): rho, rpi, rlength
%Inputs: x, y
%To modify: <none>
%Returns: srcmms_mass
% srcmasstmp; % Define return value for srcmms_mass as % precision

% dudx; 	% Temp variable: u velocity gradient in x direction
% dvdy;  % Temp variable: v velocity gradient in y direction

% This function returns the MMS mass source term

global rpi rlength
global phix phiy phixy apx apy apxy

dudx = phix(2)*apx(2)*rpi/rlength*cos(apx(2)*rpi*x/rlength)  ...
    + phixy(2)*apxy(2)*rpi*y/rlength/rlength  ...
    * cos(apxy(2)*rpi*x*y/rlength/rlength);

dvdy = -phiy(3)*apy(3)*rpi/rlength*sin(apy(3)*rpi*y/rlength)  ...
    - phixy(3)*apxy(3)*rpi*x/rlength/rlength  ...
    * sin(apxy(3)*rpi*x*y/rlength/rlength);

srcmasstmp = rho*dudx + rho*dvdy;
end
%************************************************************************
function [srcxmtmtmp] = srcmms_xmtm(x, y, rho)
%
%Uses global variable(s): rho, rpi, rmu, rlength
%Inputs: x, y
%To modify: <none>
%Returns: srcmms_xmtm

% srcxmtmtmp; % Define return value for srcmms_xmtm as % precision

% dudx; 	% Temp variable: u velocity gradient in x direction
% dudy;  % Temp variable: u velocity gradient in y direction
% termx;        % Temp variable
% termy;        % Temp variable
% termxy;       % Temp variable
% uvel;         % Temp variable: u velocity
% vvel;         % Temp variable: v velocity
% dpdx;         % Temp variable: pressure gradient in x direction
% d2udx2;       % Temp variable: 2nd derivative of u velocity in x direction
% d2udy2;       % Temp variable: 2nd derivative of u velocity in y direction

%This function returns the MMS x-momentum source term

global rpi rmu rlength
global phi0 phix phiy phixy apx apy apxy

termx = phix(2)*sin(apx(2)*rpi*x/rlength);
termy = phiy(2)*cos(apy(2)*rpi*y/rlength);
termxy = phixy(2)*sin(apxy(2)*rpi*x*y/rlength/rlength);
uvel = phi0(2) + termx + termy + termxy;

termx = phix(3)*cos(apx(3)*rpi*x/rlength);
termy = phiy(3)*cos(apy(3)*rpi*y/rlength);
termxy = phixy(3)*cos(apxy(3)*rpi*x*y/rlength/rlength);
vvel = phi0(3) + termx + termy + termxy;

dudx = phix(2)*apx(2)*rpi/rlength*cos(apx(2)*rpi*x/rlength) ...
    + phixy(2)*apxy(2)*rpi*y/rlength/rlength  ...
    * cos(apxy(2)*rpi*x*y/rlength/rlength);

dudy = -phiy(2)*apy(2)*rpi/rlength*sin(apy(2)*rpi*y/rlength)  ...
    + phixy(2)*apxy(2)*rpi*x/rlength/rlength  ...
    * cos(apxy(2)*rpi*x*y/rlength/rlength);

dpdx = -phix(1)*apx(1)*rpi/rlength*sin(apx(1)*rpi*x/rlength) ...
    + phixy(1)*apxy(1)*rpi*y/rlength/rlength  ...
    * cos(apxy(1)*rpi*x*y/rlength/rlength);

d2udx2 = -phix(2)*((apx(2)*rpi/rlength).^2)  ...
    * sin(apx(2)*rpi*x/rlength)  ...
    - phixy(2)*((apxy(2)*rpi*y/rlength/rlength).^2)  ...
    * sin(apxy(2)*rpi*x*y/rlength/rlength);

d2udy2 = -phiy(2)*((apy(2)*rpi/rlength).^2)  ...
    * cos(apy(2)*rpi*y/rlength)  ...
    - phixy(2)*((apxy(2)*rpi*x/rlength/rlength).^2)  ...
    * sin(apxy(2)*rpi*x*y/rlength/rlength);

srcxmtmtmp = rho*uvel*dudx + rho*vvel*dudy + dpdx  ...
    - rmu*( d2udx2 + d2udy2 );

end
%************************************************************************
function [srcymtmtmp] = srcmms_ymtm(x, y, rho)
%
%Uses global variable(s): rho, rpi, rmu, rlength
%Inputs: x, y
%To modify: <none>
%Returns: srcmms_ymtm

% srcymtmtmp; % Define return value for srcmms_ymtm as % precision

% dvdx;         % Temp variable: v velocity gradient in x direction
% dvdy;         % Temp variable: v velocity gradient in y direction
% termx;        % Temp variable
% termy;        % Temp variable
% termxy;       % Temp variable
% uvel;         % Temp variable: u velocity
% vvel;         % Temp variable: v velocity
% dpdy;         % Temp variable: pressure gradient in y direction
% d2vdx2;       % Temp variable: 2nd derivative of v velocity in x direction
% d2vdy2;       % Temp variable: 2nd derivative of v velocity in y direction

% This function returns the MMS y-momentum source term

global rpi rmu rlength
global phi0 phix phiy phixy apx apy apxy

termx = phix(2)*sin(apx(2)*rpi*x/rlength);
termy = phiy(2)*cos(apy(2)*rpi*y/rlength);
termxy = phixy(2)*sin(apxy(2)*rpi*x*y/rlength/rlength);
uvel = phi0(2) + termx + termy + termxy;

termx = phix(3)*cos(apx(3)*rpi*x/rlength);
termy = phiy(3)*cos(apy(3)*rpi*y/rlength);
termxy = phixy(3)*cos(apxy(3)*rpi*x*y/rlength/rlength);
vvel = phi0(3) + termx + termy + termxy;

dvdx = -phix(3)*apx(3)*rpi/rlength*sin(apx(3)*rpi*x/rlength)  ...
    - phixy(3)*apxy(3)*rpi*y/rlength/rlength  ...
    * sin(apxy(3)*rpi*x*y/rlength/rlength);

dvdy = -phiy(3)*apy(3)*rpi/rlength*sin(apy(3)*rpi*y/rlength)  ...
    - phixy(3)*apxy(3)*rpi*x/rlength/rlength  ...
    * sin(apxy(3)*rpi*x*y/rlength/rlength);

dpdy = phiy(1)*apy(1)*rpi/rlength*cos(apy(1)*rpi*y/rlength)  ...
    + phixy(1)*apxy(1)*rpi*x/rlength/rlength  ...
    * cos(apxy(1)*rpi*x*y/rlength/rlength);

d2vdx2 = -phix(3)*((apx(3)*rpi/rlength).^2)  ...
    * cos(apx(3)*rpi*x/rlength)  ...
    - phixy(3)*((apxy(3)*rpi*y/rlength/rlength).^2)  ...
    * cos(apxy(3)*rpi*x*y/rlength/rlength);

d2vdy2 = -phiy(3)*((apy(3)*rpi/rlength).^2)  ...
    * cos(apy(3)*rpi*y/rlength)  ...
    - phixy(3)*((apxy(3)*rpi*x/rlength/rlength).^2)  ...
    * cos(apxy(3)*rpi*x*y/rlength/rlength);

srcymtmtmp = rho*uvel*dvdx + rho*vvel*dvdy + dpdy  ...
    - rmu*( d2vdx2 + d2vdy2 );

end
%************************************************************************
function [dtmin] = compute_time_step(dtmin, imax, jmax, rho)
%
%Uses global variable(s): one, two, four, half, fourth
%Uses global variable(s): vel2ref, rmu, rho, dx, dy, cfl, rkappa, imax, jmax
%Uses: u
%To Modify: dt, dtmin

% i                        % i index (x direction)
% j                        % j index (y direction)

% dtvisc       % Viscous time step stability criteria (constant over domain)
% uvel2        % Local velocity squared
% beta2        % Beta squared paramete for time derivative preconditioning
% lambda_x     % Max absolute value eigenvalue in (x,t)
% lambda_y     % Max absolute value eigenvalue in (y,t)
% lambda_max   % Max absolute value eigenvalue (used in convective time step computation)
% dtconv       % Local convective time step restriction

global four half fourth
global vel2ref rmu dx dy cfl rkappa
global u dt

dtmin = +1.0e99;
dtvisc = fourth*dx*dy*rho/rmu;
%$$$$$$   fprintf('dtvisc: %f \n', dtvisc);

j=2:jmax-1;
i=2:imax-1;
uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
%$$$$$$ beta2 = rkappa*uvel2;
beta2 = max( uvel2, rkappa.*vel2ref );
lambda_x = half*( abs(u(i,j,2)) + sqrt((u(i,j,2).^2)+four.*beta2) );
lambda_y = half*( abs(u(i,j,3)) + sqrt((u(i,j,3).^2)+four.*beta2) );
lambda_max = max(lambda_x,lambda_y);
%$$$$$$     lambda_max = half*( sqrt(uvel2) + sqrt(uvel2 + four*beta2) );
dtconv = min(dx,dy) ./ ( lambda_max );
dt(i,j) = cfl*min(dtconv,dtvisc);
%$$$$$$     fprintf('dtconv, dt: %f %f\n', dtconv, dt(i,j));
%$$$$$$     dt(i,j) = fourth/( one/dtconv + one/dtvisc);
%$$$$$$     fprintf('n,i,j,dt: %d %d %d %f',n, i, j, dt(i,j));
dtmin = min(min(min(dtmin,dt(i,j))));
%dtmin = (min(dtmin,dt(i,j)));   % use in for loop version

%   UNCOMMENT BELOW FOR Global Time Stepping
%$$$$$$   for(j=1; j<jmax-1; j++)
%$$$$$$   {
%$$$$$$     for(i=1; i<imax-1; i++)
%$$$$$$     {
%$$$$$$       dt(i,j) = *dtmin;
%$$$$$$     }
%$$$$$$   }
end
%************************************************************************
function Compute_Artificial_Viscosity(imax, jmax, rho)
%
%Uses global variable(s): zero, one, two, four, six, half, fourth
%Uses global variable(s): imax, jmax, lim, rho, dx, dy, Cx, Cy, Cx2, Cy2, fsmall, vel2ref, rkappa
%Uses: u
%To Modify: artviscx, artviscy

% i                        % i index (x direction)
% j                        % j index (y direction)

% uvel2        % Local velocity squared
% beta2        % Beta squared paramete for time derivative preconditioning
% lambda_x     % Max absolute value e-value in (x,t)
% lambda_y     % Max absolute value e-value in (y,t)
% d4pdx4       % 4th derivative of pressure w.r.t. x
% d4pdy4       % 4th derivative of pressure w.r.t. y
% d2pdx2       % 2nd derivative of pressure w.r.t. x
% d2pdy2       % 2nd derivative of pressure w.r.t. y
% pfunct1      % Temporary variable for 2nd derivative damping
% pfunct2      % Temporary variable for 2nd derivative damping

global two four six half
global lim dx dy Cx Cy Cx2 Cy2 fsmall vel2ref rkappa
global u
global artviscx artviscy

% Fourth Derivative Damping
j = 3:jmax-2;
i = 3:imax-2;
d4pdx4 = ( u(i+2,j,1) - four*u(i+1,j,1) + six*u(i,j,1) ...
    - four*u(i-1,j,1) + u(i-2,j,1) )./((dx.^4));
d4pdy4 = ( u(i,j+2,1) - four*u(i,j+1,1) + six*u(i,j,1) ...
    - four*u(i,j-1,1) + u(i,j-2,1) )./((dy.^4));
uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
beta2 = max( uvel2,rkappa*vel2ref );
lambda_x = half*( abs(u(i,j,2)) + sqrt((u(i,j,2).^2)+four*beta2) );
lambda_y = half*( abs(u(i,j,3)) + sqrt((u(i,j,3).^2)+four*beta2) );
%NOTE: For Jameson damping, see also Swanson, Radespiel, and Turkel AIAA Paper (1997)
artviscx(i,j) = -lambda_x.*Cx.*(dx.^3)./beta2.*d4pdx4;   %Jameson-Turkel Damping (1981)
artviscy(i,j) = -lambda_y.*Cy.*(dy.^3)./beta2.*d4pdy4;   %Jameson-Turkel Damping (1981)
%$$$$$$     artviscx(i,j) = fourth*Cx*dx*dx*d4pdx4;  %Sotiropoulos and Abdallah (1991)
%$$$$$$     artviscy(i,j) = fourth*Cy*dy*dy*d4pdy4;  %Sotiropoulos and Abdallah (1991)




%  Side Walls
j = 3:jmax - 2;
i = 2; % Left side
d4pdy4 = ( u(i,j+2,1) - four*u(i,j+1,1) + six*u(i,j,1) ...
    - four*u(i,j-1,1) + u(i,j-2,1) )./((dy.^4));

uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
beta2 = max( uvel2,rkappa*vel2ref );
lambda_x = half*( abs(u(i,j,2)) + sqrt((u(i,j,2).^2)+four*beta2) );
lambda_y = half*( abs(u(i,j,3)) + sqrt((u(i,j,3).^2)+four*beta2) );

artviscy(i,j) = -lambda_y*Cy*(dy.^3)./beta2.*d4pdy4;  % Jameson-Turkel Damping (1981)
artviscx(i,j) = artviscx(i+1,j);

i = imax-1; % Right side
d4pdy4 = ( u(i,j+2,1) - four*u(i,j+1,1) + six*u(i,j,1) ...
    - four*u(i,j-1,1) + u(i,j-2,1) )./((dy.^4));

uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
beta2 = max( uvel2,rkappa*vel2ref );
lambda_x = half*( abs(u(i,j,2)) + sqrt((u(i,j,2).^2)+four*beta2) );
lambda_y = half*( abs(u(i,j,3)) + sqrt((u(i,j,3).^2)+four*beta2) );

artviscy(i,j) = -lambda_y*Cy*(dy.^3)./beta2.*d4pdy4;  % Jameson-Turkel Damping (1981)
artviscx(i,j) = artviscx(i-1,j);


% Top and Bottom Walls
i = 3:imax - 2;
j = 2; % Bottom
d4pdx4 = ( u(i+2,j,1) - four*u(i+1,j,1) + six*u(i,j,1) ...
    - four*u(i-1,j,1) + u(i-2,j,1) )./((dx.^4));

uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
beta2 = max( uvel2,rkappa*vel2ref );
lambda_x = half*( abs(u(i,j,2)) + sqrt((u(i,j,2).^2)+four*beta2) );
lambda_y = half*( abs(u(i,j,3)) + sqrt((u(i,j,3).^2)+four*beta2) );

artviscx(i,j) = -lambda_x*Cx*(dx.^3)./beta2.*d4pdx4;  %Jameson-Turkel Damping (1981)
artviscy(i,j) = artviscy(i,j+1);

j = jmax-1; % Top
d4pdx4 = ( u(i+2,j,1) - four*u(i+1,j,1) + six*u(i,j,1) ...
    - four*u(i-1,j,1) + u(i-2,j,1) )./((dx.^4));

uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
beta2 = max( uvel2,rkappa*vel2ref );
lambda_x = half*( abs(u(i,j,2)) + sqrt((u(i,j,2).^2)+four*beta2) );
lambda_y = half*( abs(u(i,j,3)) + sqrt((u(i,j,3).^2)+four*beta2) );

artviscx(i,j) = -lambda_x*Cx*(dx.^3)./beta2.*d4pdx4;  %Jameson-Turkel Damping (1981)
artviscy(i,j) = artviscy(i,j-1);


% Corners
i = 2;
j = 2;

artviscx(i,j) = half*(artviscx(i+1,j)+artviscx(i,j+1));
artviscy(i,j) = half*(artviscy(i+1,j)+artviscy(i,j+1));

i = 2;
j = jmax-1;

artviscx(i,j) = half*(artviscx(i+1,j)+artviscx(i,j-1));
artviscy(i,j) = half*(artviscy(i+1,j)+artviscy(i,j-1));

i = imax-1;
j = 2;

artviscx(i,j) = half*(artviscx(i-1,j)+artviscx(i,j+1));
artviscy(i,j) = half*(artviscy(i-1,j)+artviscy(i,j+1));


i = imax-1;
j = jmax-1;
artviscx(i,j) = half*(artviscx(i-1,j)+artviscx(i,j-1));
artviscy(i,j) = half*(artviscy(i-1,j)+artviscy(i,j-1));

% Second Derivative Damping
if ((Cy2~=0.0)&&(Cx2~=0.0))
    for j=2:jmax-1
        for i=2:imax-1
            pfunct1 = abs( u(i+1,j,lim) - two*u(i,j,lim) + u(i-1,j,lim) );
            pfunct2 = abs( u(i+1,j,lim) + two*u(i,j,lim) + u(i-1,j,lim) );
            
            d2pdx2 = ( u(i+1,j,lim) - two*u(i,j,lim) + u(i-1,j,lim) )/(dx*dx);
            artviscx(i,j) = artviscx(i,j) ...
                + rho*Cx2*dx*d2pdx2*pfunct1/(pfunct2+fsmall);  %2nd form
            
            pfunct1 = abs( u(i,j+1,lim) - two*u(i,j,lim) + u(i,j-1,lim) );
            pfunct2 = abs( u(i,j+1,lim) + two*u(i,j,lim) + u(i,j-1,lim) );
            
            d2pdy2 = (u(i,j+1,lim) - two*u(i,j,lim) + u(i,j-1,lim) )/(dy*dy);
            artviscy(i,j) = artviscy(i,j)  ...
                + rho*Cy2*dy*d2pdy2*pfunct1/(pfunct2+fsmall);   %2nd form
        end
    end
end

end
%************************************************************************
function SGS_forward_sweep(imax, jmax, rho)
%
%Uses global variable(s): two, three, six, half
%Uses global variable(s): imax, imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, ...
%                      xmax, xmin, ymax, ymin, rmu, vel2ref
%Uses: artviscx, artviscy, dt, s
%To Modify: u

% i                        % i index (x direction)
% j                        % j index (y direction)

% dpdx         % First derivative of pressure w.r.t. x
% dudx         % First derivative of x velocity w.r.t. x
% dvdx         % First derivative of y velocity w.r.t. x
% dpdy         % First derivative of pressure w.r.t. y
% dudy         % First derivative of x velocity w.r.t. y
% dvdy         % First derivative of y velocity w.r.t. y
% d2udx2       % Second derivative of x velocity w.r.t. x
% d2vdx2       % Second derivative of y velocity w.r.t. x
% d2udy2       % Second derivative of x velocity w.r.t. y
% d2vdy2       % Second derivative of y velocity w.r.t. y
% beta2        % Beta squared parameter for time derivative preconditioning
% uvel2        % Velocity squared

global two half
global rhoinv dx dy rkappa rmu vel2ref
global artviscx artviscy dt s u
%------------------------------------------------------------------------------------------------------------------%
global rhs
%------------------------------------------------------------------------------------------------------------------%

% Symmetric Gauss-Siedel: Forward Sweep
for i = 2:imax-1
    for j = 2:jmax-1
        dpdx = half.*(u(i+1,j,1) - u(i-1,j,1))./dx;
        dudx = half.*(u(i+1,j,2) - u(i-1,j,2))./dx;
        dvdx = half.*(u(i+1,j,3) - u(i-1,j,3))./dx;
        dpdy = half.*(u(i,j+1,1) - u(i,j-1,1))./dy;
        
        dudy = half*(u(i,j+1,2) - u(i,j-1,2))./dy;
        dvdy = half*(u(i,j+1,3) - u(i,j-1,3))./dy;
        d2udx2 = (u(i+1,j,2) - two*u(i,j,2) + u(i-1,j,2))./(dx.*dx);
        d2vdx2 = (u(i+1,j,3) - two*u(i,j,3) + u(i-1,j,3))./(dx.*dx);
        d2udy2 = (u(i,j+1,2) - two*u(i,j,2) + u(i,j-1,2))./(dy.*dy);
        d2vdy2 = (u(i,j+1,3) - two*u(i,j,3) + u(i,j-1,3))./(dy.*dy);
        uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
        beta2 = max( uvel2,rkappa.*vel2ref );
        %----------------------------------------------------------------------------------------------------------%
        rhs(i,j,1) = beta2.*( rho.*dudx + rho.*dvdy  ...
            - artviscx(i,j) - artviscy(i,j)  ...
            - s(i,j,1) );
        rhs(i,j,2) = ( u(i,j,2).*dudx + u(i,j,3).*dudy ...
            + rhoinv.*dpdx  ...
            - rhoinv.*rmu.*d2udx2 - rhoinv.*rmu.*d2udy2  ...
            - rhoinv.*s(i,j,2) );
        rhs(i,j,3) = ( u(i,j,2).*dvdx + u(i,j,3).*dvdy + rhoinv.*dpdy  ...
            - rhoinv.*rmu.*d2vdx2 - rhoinv.*rmu.*d2vdy2 ...
            - rhoinv.*s(i,j,3) );
        
        u(i,j,1) = u(i,j,1)  ...
            - dt(i,j).*rhs(i,j,1);
        u(i,j,2) = u(i,j,2)  ...
            - dt(i,j).*rhs(i,j,2);
        u(i,j,3) = u(i,j,3) ...
            - dt(i,j).*rhs(i,j,3);
        %----------------------------------------------------------------------------------------------------------%
    end
end


end
%************************************************************************
function SGS_backward_sweep(imax, jmax, rho)
%
%Uses global variable(s): two, three, six, half
%Uses global variable(s): imax, imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, ...
%                      xmax, xmin, ymax, ymin, rmu, vel2ref
%Uses: artviscx, artviscy, dt, s
%To Modify: u

% i                        % i index (x direction)
% j                        % j index (y direction)

% dpdx         % First derivative of pressure w.r.t. x
% dudx         % First derivative of x velocity w.r.t. x
% dvdx         % First derivative of y velocity w.r.t. x
% dpdy         % First derivative of pressure w.r.t. y
% dudy         % First derivative of x velocity w.r.t. y
% dvdy         % First derivative of y velocity w.r.t. y
% d2udx2       % Second derivative of x velocity w.r.t. x
% d2vdx2       % Second derivative of y velocity w.r.t. x
% d2udy2       % Second derivative of x velocity w.r.t. y
% d2vdy2       % Second derivative of y velocity w.r.t. y
% beta2        % Beta squared parameter for time derivative preconditioning
% uvel2        % Velocity squared

global two half
global rhoinv dx dy rkappa rmu vel2ref
global artviscx artviscy dt s u
%------------------------------------------------------------------------------------------------------------------%
global rhs
%------------------------------------------------------------------------------------------------------------------%

% Symmetric Gauss-Siedel: Backward Sweep
for j=jmax-1 : -1 : 2
    for i=imax-1 : -1 : 2
        dpdx = half*(u(i+1,j,1) - u(i-1,j,1))/dx;
        dudx = half*(u(i+1,j,2) - u(i-1,j,2))/dx;
        dvdx = half*(u(i+1,j,3) - u(i-1,j,3))/dx;
        dpdy = half*(u(i,j+1,1) - u(i,j-1,1))/dy;
        
        dudy = half*(u(i,j+1,2) - u(i,j-1,2))/dy;
        dvdy = half*(u(i,j+1,3) - u(i,j-1,3))/dy;
        d2udx2 = (u(i+1,j,2) - two*u(i,j,2) + u(i-1,j,2))/(dx*dx);
        d2vdx2 = (u(i+1,j,3) - two*u(i,j,3) + u(i-1,j,3))/(dx*dx);
        d2udy2 = (u(i,j+1,2) - two*u(i,j,2) + u(i,j-1,2))/(dy*dy);
        d2vdy2 = (u(i,j+1,3) - two*u(i,j,3) + u(i,j-1,3))/(dy*dy);
        uvel2 = (u(i,j,2).^2) + (u(i,j,3).^2);
        beta2 = max( uvel2,rkappa*vel2ref );
        
        %----------------------------------------------------------------------------------------------------------%
        rhs(i,j,1) = beta2.*( rho.*dudx + rho.*dvdy  ...
            - artviscx(i,j) - artviscy(i,j)  ...
            - s(i,j,1) );
        rhs(i,j,2) = ( u(i,j,2).*dudx + u(i,j,3).*dudy ...
            + rhoinv.*dpdx  ...
            - rhoinv.*rmu.*d2udx2 - rhoinv.*rmu.*d2udy2  ...
            - rhoinv.*s(i,j,2) );
        rhs(i,j,3) = ( u(i,j,2).*dvdx + u(i,j,3).*dvdy + rhoinv.*dpdy  ...
            - rhoinv.*rmu.*d2vdx2 - rhoinv.*rmu.*d2vdy2 ...
            - rhoinv.*s(i,j,3) );
        
        u(i,j,1) = u(i,j,1)  ...
            - dt(i,j).*rhs(i,j,1);
        u(i,j,2) = u(i,j,2)  ...
            - dt(i,j).*rhs(i,j,2);
        u(i,j,3) = u(i,j,3) ...
            - dt(i,j).*rhs(i,j,3);
        %----------------------------------------------------------------------------------------------------------%
    end
end
end
%************************************************************************
function point_Jacobi(imax, jmax, rho)
%
%Uses global variable(s): two, three, six, half
%Uses global variable(s): imax, imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, ...
%                      xmax, xmin, ymax, ymin, rmu, vel2ref
%Uses: uold, artviscx, artviscy, dt, s
%To Modify: u


% i                        % i index (x direction)
% j                        % j index (y direction)

% dpdx         % First derivative of pressure w.r.t. x
% dudx         % First derivative of x velocity w.r.t. x
% dvdx         % First derivative of y velocity w.r.t. x
% dpdy         % First derivative of pressure w.r.t. y
% dudy         % First derivative of x velocity w.r.t. y
% dvdy         % First derivative of y velocity w.r.t. y
% d2udx2       % Second derivative of x velocity w.r.t. x
% d2vdx2       % Second derivative of y velocity w.r.t. x
% d2udy2       % Second derivative of x velocity w.r.t. y
% d2vdy2       % Second derivative of y velocity w.r.t. y
% beta2        % Beta squared parameter for time derivative preconditioning
% uvel2        % Velocity squared
global two half
global rhoinv dx dy rkappa rmu vel2ref
global u uold artviscx artviscy dt s
%------------------------------------------------------------------------------------------------------------------%
global rhs
%------------------------------------------------------------------------------------------------------------------%

% Point Jacobi method
j=2:jmax-1;
i=2:imax-1;
dpdx = half*(uold(i+1,j,1) - uold(i-1,j,1))./dx;
dudx = half*(uold(i+1,j,2) - uold(i-1,j,2))./dx;
dvdx = half*(uold(i+1,j,3) - uold(i-1,j,3))./dx;
dpdy = half*(uold(i,j+1,1) - uold(i,j-1,1))./dy;

dudy = half*(uold(i,j+1,2) - uold(i,j-1,2))./dy;
dvdy = half*(uold(i,j+1,3) - uold(i,j-1,3))./dy;
d2udx2 = (uold(i+1,j,2) - two*uold(i,j,2) + uold(i-1,j,2))./(dx.*dx);
d2vdx2 = (uold(i+1,j,3) - two*uold(i,j,3) + uold(i-1,j,3))./(dx.*dx);
d2udy2 = (uold(i,j+1,2) - two*uold(i,j,2) + uold(i,j-1,2))./(dy.*dy);
d2vdy2 = (uold(i,j+1,3) - two*uold(i,j,3) + uold(i,j-1,3))./(dy.*dy);
uvel2 = (uold(i,j,2).^2) + (uold(i,j,3).^2);
beta2 = max( uvel2,rkappa.*vel2ref );

%----------------------------------------------------------------------------------------------------------%
rhs(i,j,1) = beta2.*( rho.*dudx + rho.*dvdy  ...
    - artviscx(i,j) - artviscy(i,j)  ...
    - s(i,j,1) );
rhs(i,j,2) = ( u(i,j,2).*dudx + u(i,j,3).*dudy ...
    + rhoinv.*dpdx  ...
    - rhoinv.*rmu.*d2udx2 - rhoinv.*rmu.*d2udy2  ...
    - rhoinv.*s(i,j,2) );
rhs(i,j,3) = ( u(i,j,2).*dvdx + u(i,j,3).*dvdy + rhoinv.*dpdy  ...
    - rhoinv.*rmu.*d2vdx2 - rhoinv.*rmu.*d2vdy2 ...
    - rhoinv.*s(i,j,3) );
        
u(i,j,1) = uold(i,j,1)  ...
    - dt(i,j).*rhs(i,j,1);
u(i,j,2) = uold(i,j,2)  ...
    - dt(i,j).*rhs(i,j,2);
u(i,j,3) = uold(i,j,3) ...
    - dt(i,j).*rhs(i,j,3);
%----------------------------------------------------------------------------------------------------------%

end
%************************************************************************
function pressure_rescaling(imax, jmax, imms)
%
%Uses global variable(s): imax, jmax, imms, xmax, xmin, ymax, ymin, rlength, pinf
%To Modify: u

% i                        % i index (x direction)
% j                        % j index (y direction)

% iref                      % i index location of pressure rescaling point
% jref                      % j index location of pressure rescaling point

% x        % Temporary variable for x location
% y        % Temporary variable for y location
% deltap   % delta_pressure for rescaling all values

global xmax xmin ymax ymin pinf
global u

iref = (imax-1)/2+1;     % Set reference pressure to center of cavity
jref = (jmax-1)/2+1;
if (imms==1)
    x = (xmax - xmin)*(iref-1)/(imax - 1);
    y = (ymax - ymin)*(jref-1)/(jmax - 1);
    deltap = u(iref,jref,1) - umms(x,y,1); % Constant in MMS
else
    deltap = u(iref,jref,1) - pinf; % Reference pressure
end

j=1:jmax;
i=1:imax;
u(i,j,1) = u(i,j,1) - deltap;


end
%************************************************************************
function [res, resinit, conv] = check_iterative_convergence...
    (n, res, resinit, ninit, rtime, dtmin, imax, jmax)
%
%Uses global variable(s): zero
%Uses global variable(s): imax, jmax, neq, fsmall
%Uses: n, u, uold, dt, res, resinit, ninit, rtime, dtmin
%To modify: conv

% i                        % i index (x direction)
% j                        % j index (y direction)
% k                        % k index (# of equations)

global zero
global neq fsmall
global u uold dt fp1

% Compute iterative residuals to monitor iterative convergence
k =1:neq;
res(k) = zero;
j = 2:jmax-1;
i = 2:imax-1;
res(1) = res(1) + sum(sum( ((u(i,j,1) - uold(i,j,1))./dt(i,j) ).^2 ));
res(2) = res(2) + sum(sum( ((u(i,j,2) - uold(i,j,2))./dt(i,j) ).^2 ));
res(3) = res(3) + sum(sum( ((u(i,j,3) - uold(i,j,3))./dt(i,j) ).^2 ));
% Should be OK for explicit schemes, but need to check
res(k) = sqrt(res(k)./( (imax-2).*(jmax-2) ));
if (n<6)
    resinit(k) = max(res(k),fsmall);
    %$$$$$$       printf("resinit, k: %f %d\n",resinit(k), k);
end
res(k) = res(k)./resinit(k);


conv = max(res);

% Write iterative residuals every 10 iterations
if ( (mod(n,10)==0)||(n==ninit) )
    fprintf(fp1, '%d %e %e %e %e\n',n, rtime, res(1), res(2), res(3) );
    fprintf('%d   %e   %e   %e   %e   %e\n',n, rtime, dtmin, res(1), res(2), res(3) );
    % Maybe a need to format this better
end

% Write header for iterative residuals every 200 iterations
if ( (mod(n,200)==0)||(n==ninit) )
    fprintf('Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum\n');
end

end

%************************************************************************
function Discretization_Error_Norms(rL1norm, rL2norm, rLinfnorm, imax, jmax, imms)
%
%Uses global variable(s): zero
%Uses global variable(s): imax, jmax, neq, imms, xmax, xmin, ymax, ymin, rlength
%Uses: u
%To modify: rL1norm, rL2norm, rLinfnorm


% i                        % i index (x direction)
% j                        % j index (y direction)
% k                        % k index (# of equations)

% x        % Temporary variable for x location
% y        % Temporary variable for y location
% DE   	% Discretization error (absolute value)

global zero neq xmax xmin ymax ymin u

if imms==1
    for k=1:neq
        rL1norm(k) = zero;
        rL2norm(k) = zero;
        rLinfnorm(k) = zero;
        for j=1:jmax
            for i=1:imax
                x = (xmax - xmin)*(i-1)/(imax - 1);
                y = (ymax - ymin)*(j-1)/(jmax - 1);
                DE = abs(u(i,j,k) - umms(x,y,k));
                rL1norm(k) = rL1norm(k) + DE;
                rL2norm(k) = rL2norm(k) + DE*DE;
                rLinfnorm(k) = max(rLinfnorm(k),DE);
            end
        end
        rL1norm(k) = rL1norm(k)/(imax*jmax);
        rL2norm(k) = sqrt(rL2norm(k)/(imax*jmax) );
    end
    
    % Write DE norms output file
    fp5 = fopen('./DE_norms.dat','w');
    fprintf('\n');
    fprintf('#################################');
    fprintf('\n');
    fprintf('DE Norms for grid: %d x %d\n', imax, jmax);
    fprintf(fp5,'DE Norms for grid: %d x %d\n', imax, jmax);
    fprintf('Variable      L1norm        L2norm        Linfnorm\n');
    fprintf(fp5,'Variable      L1norm        L2norm        Linfnorm\n');
    for k=1:neq
        fprintf('%d  %e  %e  %e\n', k, rL1norm(k), rL2norm(k), rLinfnorm(k));
        fprintf(fp5, '%d  %e  %e  %e\n', k, rL1norm(k), rL2norm(k), rLinfnorm(k));
    end
    fclose(fp5);
end

end
