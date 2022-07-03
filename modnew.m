function [a,iter]=modnew(tol,ak,rk,Vr,Q,rstr,...
    uinf, rho, imax, jmax, imms, isgs, uref,...
    loc_p, loc_u, loc_v)

%Function to arrive at the converged (steady) generalized co-ordinates
%using modified Gauss-Newton 
%Inputs: Steady state residual tolerance, tol 
%       Small increment, e
%       Generalized co-ordinates ak (neq*rx1)
%       Reference state, ubar (neq*Nx1)
%       Residuals, rk (neq*rx1)
%       POD Basis, Vr (neq*Nxneq*r)
%       RHS POD Basis, Q (neq*rxneq*N)
%       Restart condition, 0 (no restart), 1 (use restart data)
% Rest of the inputs as required by the cavity_solver_iter function
%Outputs: Converged generalized co-ordinates: a(rx1)

%% Set plot stuff:
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultLineColor',[1,1,1])
set(0,'DefaultLineMarkerSize',15)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultFigureColor',[1,1,1])
set(0,'DefaultTextFontSize',18)
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultTextFontName','Times-Roman')
set(0,'DefaultAxesFontName','Times-Roman')

%%
% Maximum number of iterations
maxiter = 50;
% Output write interval
interval = 5;
% No. of equations (p,u,v)
neq = 3;
% Initialize the residual norm vector
d0 = zeros(neq,1);

% Full order dimension
N = imax*jmax; % size(Vr,1)/neq;
% Reduced order dimension
r = length(ak)/neq;

switch rstr
    case 0
        
        % Initial residual norms for individual scaling
        % Scale the residuals with their respective initial norms if the
        % initial are not zero
        for i = 1 : neq
            if norm(rk((i-1)*r+1:i*r),2)>1e-12
                rk((i-1)*r+1:i*r) = rk((i-1)*r+1:i*r)/norm(rk((i-1)*r+1:i*r),2);
                d0(i) = norm(rk((i-1)*r+1:i*r),2);
            end
        end
        % Initial residual norms vector
        dk = d0;
        % Record Convergence history
        Dk = dk;
        % Set iteration counter
        iter = 1;
    case 1
        % Load saved variables
        load('hrom_restart','iter','ak','rk','Dk','Vr','Q');
        % Initial residual norms for individual scaling
        d0 = Dk(:,1);
        %Residual norms vector
        dk = Dk(:,end);
end

%Initialize total timer
totalitert = 0; %min

while all(dk>tol)||iter<maxiter
    %Iteration timer
    tstart = tic;
    %Iteration counter
    iter = iter+1;
    
    % Modified Gauss-Newton method 
    
    % Obtain mu(x) = f(x)/f'(x); mu(x) = J\f(x)
    
    mu = getmu(neq,r,ak,rk,Vr,Q,...
    uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0);
    
    % Determine dak = -mu(x)/mu'(x)
    Jmu = fdjacmodnew(mu,ak,rk,Vr,Q,...
          uinf, rho, imax, jmax, imms, isgs, uref,...
          loc_p, loc_u, loc_v, d0);
        
    % Initiate equation index
    k = 1;
    while k <= 3 && norm(rk((k-1)*r+1:k*r),2)>tol
        if rank(Jmu((k-1)*r+1:k*r,(k-1)*r+1:k*r)) ~= length(ak)
            dak = ...
            -pinv(Jmu((k-1)*r+1:k*r,(k-1)*r+1:k*r))*mu((k-1)*r+1:k*r);
        else
            dak = ...
            -Jmu((k-1)*r+1:k*r,(k-1)*r+1:k*r)\mu((k-1)*r+1:k*r);
        end
        % Update reduced co-ordinates
        ak((k-1)*r+1:k*r) = ak((k-1)*r+1:k*r)+dak;
        % Update k
        k = k+1;
    end
    
    % Get RHS
    rk = getrhs(uref,Vr,ak,Q,uinf, rho, imax, jmax, imms, isgs,...
    loc_p, loc_u, loc_v);
    
    %Initial residual norms for individual scaling
    %Scale the residuals with their respective initial norms
    for i = 1 : neq
        if norm(rk((i-1)*r+1:i*r),2)>1e-12
                rk((i-1)*r+1:i*r) = rk((i-1)*r+1:i*r)/d0(i);
                dk(i) = norm(rk((i-1)*r+1:i*r),2);
        end
    end
    
    %Record Convergence history
    Dk = [Dk,dk];
    save('hrom_conv.mat','Dk');
    
    if iter<maxiter && mod(iter,interval)==0
        save('hrom_restart','iter','ak','rk','Dk','Vr','Q');
    elseif iter==maxiter
        fprintf('Max. number of iterations reached before convergence\n');
        break;
    end
       
    %Iteration timer
    oneitert=toc(tstart);
    oneitertm=fix(oneitert/60); oneiterts=mod(oneitert,60); %min, sec
    totalitert=totalitert+oneitert;
    totalitertm=fix(totalitert/60);totaliterts=mod(totalitert,60); %min, sec
    %Iteration display
    fprintf('Iteration: %d\t time: %d m %d s\t total time: %d m % s\n',...
        iter,oneitertm,oneiterts,totalitertm,totaliterts);
    fprintf('p_res: %.8f\t, u_res: %.8f\t, v_res: %.8f\n',...
        dk(1),dk(2),dk(3));
end

a = ak;


% Plot convergence history
figure(1)
semilogy(iter,Dk(1,:),iter,Dk(2,:),iter,Dk(3,:))
xlabel('Iterations')
ylabel('$||R||_2$')
xlim([1,50])
legend('p','u','v')
grid on
    

end
