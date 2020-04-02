function [Vf] = polyfilt(V,sten,inCyl)
% Function that performs a Savitzky-Golay spatial filter of 2nd order. The
% border points are not filtered.
%
% INPUT:
%   V     - Unfiltered velocity field, one component. Size (nx,ny)
%   sten  - Stencil size used in the spatial filter [pixels]
%   inCyl - Binary matrix that says wether the grid point is inside the cylinder or not
%
% OUTPUT:
%   Vf    - Filtered velocity field, one component. Size (nx,ny)

% Calculate matrix coefficients
cont = 0;                                   % Counter
M = zeros(2*length(-sten:sten),6);          % Pre-allocation for memory
for i = -sten:sten                          % x coordinate
    for j = -sten:sten                      % y coordinate
        cont = cont+1;                      % Update counter
        M(cont,:) = [1 j i j*j j*i i*i];    % Coefficients
    end
end
Mat=(M'*M)\(M');

% Pre-allocation for memory
[A,B] = size(V);
Vf    = V;

% Filtering process
for i=(sten+1):A-sten           % For each row
    for j=(sten+1):B-sten       % For each column
        if ~inCyl(i,j)          % If central point is outside the cylinder
            c1 = i-sten:i+sten;
            c2 = j-sten:j+sten;
            dum=V(i-sten:i+sten,j-sten:j+sten);
            
            % Check if any of the stencil points are inside the cylinder
            for ii = 1:length(i-sten:i+sten)
                for jj = 1:length(j-sten:j+sten)
                    if inCyl(c1(ii),c2(jj))
                        dum(ii,jj) = 0;     % If they are, neglect them
                    end
                end
            end
            
            % Multiply window with matrix coefficients
            dum = dum(:);
            app = Mat*dum;
            
            % Extract filtered value
            Vf(i,j) = app(1);
        end
    end
end