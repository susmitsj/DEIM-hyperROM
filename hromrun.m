clc
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        DEIM-HROM Framework        %%%%%%%%%%%%
%%%%           Prepared by Susmit Joshi              %%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Code to call the hrom function and store the variables in an output file

%%

% Discretization
imax = 65;
jmax = 65;

% FOM 
imms = 1; % Manufactured solution on
isgs = 1; % Symmetric Gauss-Seidel on

% ROM Test
N = 1;

% Parameter limits
uinf_min = 0.5;   % Lid velocity (m/s)
uinf_max = 2;
rho_min = 0.5;    % Air density (kg/m^3)
rho_max = 2;

% MC Testing Parameters
% rng(0,'twister');
% seq = randn(N,1);

uinf = uinf_min + (uinf_max - uinf_min)*(randi([0,1]));

rho = rho_min + (rho_max - rho_min)*(randi([0,1]));

x = [uinf, rho];

% Simulation number

sim_no = 1;

%% QDEIM HROM

method = 'qdeim';
[rom_p, rom_u, rom_v, rom_iter] = hrom(method, x, imax, jmax, imms, isgs);

%% FOM and MS

[fom_p, fom_u, fom_v, p, u, v] = cavity_solver_mms(uinf, rho, imax, jmax, isgs, sim_no);


%%




% %% ROM
% method='no';
% [u3,iter3]=hrom(method);
% 
% outputfile = datafolder+"hrom_output3.dat";
% 
% fileID=fopen(outputfile,'w');
% fprintf(fileID,'U\n');
% for i=1:length(u1)
%     fprintf(fileID,'%.8f\n',u3(i)); 
% end
% fclose(fileID);
