%**************************************************************************
% RELAX_PSI.m
% Last edited by: pjh4 Nov 2020
%
% This function does Gauss iteration to solve the Poisson equation
%
% grid - the struct containing data on the simulation
% vort_on - is equal to 1 if vort in grid should be used (otherwise Poisson
%           eqn will degenerate to Laplace)
% over_relax - 'false' if over relaxation should not be used
% cycle - the iteration in the simulation
% iter - number of iterations needed for Poisson equation to converge
%**************************************************************************

function [grid, iter] = RELAX_PSI(grid, vort_on, over_relax, cycle)

switch over_relax
    
    case 'false'
        
        if (vort_on == 1)
            
            disp('Relaxing psi without over-relaxation, with vorticity')
            [grid, iter] = RELAX_PSI_WITH_VORT_NO_F(grid, cycle);
        else
            
            disp('Relaxing psi without over-relaxation, no vorticity')
            [grid, iter] = RELAX_PSI_NO_F(grid, cycle);
        end
        
    otherwise
        
        if (vort_on == 1)
            
            disp('Relaxing psi with over-relaxation, with vorticity')
            [grid, iter] = RELAX_PSI_WITH_VORT(grid, cycle);
        else
            
            disp('Relaxing psi with over-relaxation, no vorticity')
            [grid, iter] = RELAX_PSI_NO_VORT(grid, cycle);
        end
end

end


function [grid, k] = RELAX_PSI_NO_F(grid, cycle)

if cycle == 1
    psi = grid.psi(:,:,cycle);
else
    psi = grid.psi(:,:,cycle-1);
end

k = 1;
while 1
    %     k
    maxdiff = 0;
    for ii=2:grid.rows-1
        for j=2:grid.cols-1
            
            % store old grid.psi value
            oldValue = psi(ii,j);
            
            % update stream function values that are not in/on the object
            if (grid.key(ii,j) == 0)
                
                psi(ii,j) = 1/4*(psi(ii+1,j) + ...
                    psi(ii-1,j) + psi(ii,j+1) + ...
                    psi(ii,j-1));
                
            end
            
            % check convergence criterion
            if oldValue == 0
                temp_diff = 0;
            else
                temp_diff = abs((psi(ii,j)-oldValue)/oldValue);
            end
            
            maxdiff = max([maxdiff, temp_diff]);
            
        end
    end

    if (maxdiff < grid.epsilon)
        disp('Convergence reached for bulk computation')
        break;
    end
    
    k = k+1;
end

grid.psi(:,:,cycle) = psi;

end

function [grid, k] = RELAX_PSI_NO_VORT(grid, cycle)

if cycle == 1
    psi = grid.psi(:,:,cycle);
else
    psi = grid.psi(:,:,cycle-1);
end

k = 1;
while 1
    maxdiff = 0;
    for ii=2:grid.rows-1
        for j=2:grid.cols-1
            
            % store old grid.psi value
            oldValue = psi(ii,j);
            
            % update stream function values that are not in/on the object
            if (grid.key(ii,j) == 0)
                
                psi(ii,j) = oldValue + grid.F/4*(psi(ii+1,j) + ...
                    psi(ii-1,j) + psi(ii,j+1) + ...
                    psi(ii,j-1) - 4 * oldValue);
                
            end
            
            % check convergence criterion
            if oldValue == 0
                temp_diff = 0;
            else
                temp_diff = abs((psi(ii,j)-oldValue)/oldValue);
            end
            
            maxdiff = max([maxdiff, temp_diff]);
        end
    end

    if (maxdiff < grid.epsilon)
        disp('Convergence reached for bulk computation')
        break;
    end
    
    k = k+1;
end

grid.psi(:,:,cycle) = psi;

end

function [grid, k] = RELAX_PSI_WITH_VORT(grid, cycle)

psi = grid.psi(:,:,cycle-1);
vort = grid.vort(:,:,cycle-1);

grid_old = psi;

k = 1;
while 1

    for ii=2:grid.rows-1
        for j=2:grid.cols-1
            
            % update stream function values that are not in/on the object
            if (grid.key(ii,j) == 0)
                
                psi(ii,j) = psi(ii,j) + grid.F/4*(psi(ii+1,j) + ...
                    psi(ii-1,j) + psi(ii,j+1) + ...
                    psi(ii,j-1) + 4 * (grid.h)^2 * vort(ii,j)- 4*psi(ii,j));
                
            end
            
        end
    end
    
    temp = abs((grid_old(grid_old~=0) - psi(grid_old~=0))./grid_old(grid_old~=0));
    maxdiff = max(max(temp));
    
    grid_old = psi;
    
    if (maxdiff < grid.epsilon)
        disp('Convergence reached for bulk computation')
        break;
    end
    
    k = k+1;
    
end

grid.psi(:,:,cycle) = psi;

end


function [grid, k] = RELAX_PSI_WITH_VORT_NO_F(grid, cycle)

psi = grid.psi(:,:,cycle-1);
vort = grid.vort(:,:,cycle-1);

grid_old = psi;

k = 1;
while 1

    for ii=2:grid.rows-1
        for j=2:grid.cols-1
            
            % update stream function values that are not in/on the object
            if (grid.key(ii,j) == 0)
                
                psi(ii,j) = 1/4*(psi(ii+1,j) + ...
                    psi(ii-1,j) + psi(ii,j+1) + ...
                    psi(ii,j-1) + 4 * (grid.h)^2 * vort(ii,j));
                
            end
        end
    end
    
    temp = abs((grid_old(grid_old~=0) - psi(grid_old~=0))./grid_old(grid_old~=0));
    maxdiff = max(max(temp));
    
    grid_old = psi;
    
    if (maxdiff < grid.epsilon)
        disp('Convergence reached for bulk computation')
        break;
    end
    
    k = k+1;
    
end

grid.psi(:,:,cycle) = psi;

end

