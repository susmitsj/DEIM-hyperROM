function [Jmu] = fdjacmodnew(mu,ak,rk,Vr,Q,...
    uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0)

%Function to calculate the finite difference Jacobian of mu(x) =
%f(x)/f'(x)from the modified newton approach
%Inputs:Modified Newton's function mu(x)=f(x)/f'(x) (neq*rx1)
%       Generalized co-ordinates ak (neq*rx1)
%       Reference state, ubar (neq*Nx1)
%       Residuals, mu (neq*rx1)
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
Jmu = zeros(neq*r,neq*r); % Initialized as zeros matrix

%Small change in ak
e = 1e-2*ak;

% Parallel processing 
% p=parpool(3); % Run this before you start the code
p = gcp();
%Initiate parfeval futures
f(1:r)=parallel.FevalFuture;
%Initialize Output matrix
F=cell(r);


%Determine the Jacobian of rk w.r.t. ak
for k = 1 : neq
    for j = 1 : r
        %Get the corresponding ak to the equation
        akk = ak((k-1)*r+1:k*r);
        %Change each element by e for every j
        akk(j) = akk(j)+e((k-1)*r+j);
        %Create a (neq*rx1) from the new vector
        a = ak;
        a((k-1)*r+1:k*r) = akk;
        
        % Get mu(x)
%         muk = getmu(neq,r,a,rk,Vr,Q,...
%             uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0);    
        
        f(j)=parfeval(p,@getmu,1,neq,r,a,rk,Vr,Q,...
            uinf,rho,imax,jmax,imms,isgs,uref,loc_p,loc_u,loc_v,d0);
        
%         % Create Jacobian
%         Jmu((k-1)*r+1:k*r,(k-1)*r+j) = ...
%             (muk((k-1)*r+1:k*r)-mu((k-1)*r+1:k*r))/e((k-1)*r+j);
    end
    
    wait(f);
    cancel(p.FevalQueue.QueuedFutures);
    cancel(p.FevalQueue.RunningFutures);
    
    for j=1:r
        [completedId,value]=fetchNext(f);
        F{completedId}=value;
%             fprintf('Got result with index: %d .\n',completedId);
    end
    
    Muk=cell2mat(F);
    %Get Jacobian (neq*rxneq*r) 
    for j=1:r
        muk=Muk((j-1)*r+1:j*r);
        Jmu((k-1)*r+1:k*r,(k-1)*r+j)=(muk-mu((k-1)*r+1:k*r))/e((k-1)*r+j);
    end
    
end

end
