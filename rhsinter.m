function [loc] = rhsinter(jmax, p)

% Function to obtain the new u (imax x jmax) with non-zero entries
% only at the points given by the selection matrix P^T and the neighboring
% points as needed for the finite difference scheme in the solver function
% (i-1,i-1,i+1,i+2) and (j-2,j-1,j+1,j+2)

% Inputs
% u: Solution matrix for each equation (imax x jmax)
% p: Vector selection operator from QDEIM algorithm (mx1)
% Outputs
% unew: Solution matrix with non-zero values at necessary places
% loc: Matrix containing the locations for interpolation points (m x 2)

% Get the reduced dimension for right hand side
m = length(p);
% Obtain imax and jmax
% imax = size(u,1);
% jmax = size(u,2);

% % Initiate unew
% unew = zeros(size(u));
% Initiate loc
loc = zeros(m,2);

% Get (i, j) locations from p
for k = 1:m
    % After division by jmax, i=quotient+1, j=remainder
    loc(k,2) = mod(p(k),jmax);
    loc(k,1) = (p(k)-loc(k,2))/jmax + 1;
    
%     % Fill unew with u(i,j) based on loc and the neighboring points needed for
%     % finite differences
%     for i = 1:imax
%         for j = 1:jmax
%             if loc(k,1)-2 <= i <= loc(k,1)+2 && loc(k,2)-2 <= j <= loc(k,2)+2
%                 unew(i,j) = u(i,j);
%             end
%         end
%     end
end






