function [ phase_series, oscillation_flag ] = map_phase( cycle, y )
%MAP_PHASE Maps new time series data onto corresponding phase of limit cycle
%   Takes limit cycle and new time series data y, and matches the
%   time series data onto the limit cycle, returning a time series of the
%   phase and a flag for whether or not there are oscillations

half_length = round(3*size(y,1)/4);
tol = 1e-2;

%check if oscillation is even there by checking if derivative has a sign
%change
sign_changes=0;
for ii=half_length:size(y,1)
if sign(y(ii,1)-y(ii-1,1)) ~= sign(y(ii-1,1)-y(ii-2,1))
    display('sign changes')
    sign_changes = 1;
    break
end
end
if sign_changes==0
    disp('no oscillation')
    oscillation_flag = 0;
    phase_series = 0*y(half_length:end,1);
    return
end

%otherwise, map phase
oscillation_flag = 1;
start_flag = 0;
T_min = 1e1;
for ii=half_length:size(y,1)
   if abs(y(ii,1)) < tol && start_flag == 0 
       start_ind = ii;
       start_flag = 1;
       if y(ii,1)>y(ii-1,1)
           upward_flag = 1;
       else
           upward_flag = 0;
       end
   elseif start_flag == 1 && ii - start_ind < T_min
       continue
   elseif abs(y(ii,1)) < tol && start_flag == 1 && y(ii,1)>y(ii-1,1) && upward_flag == 1
           end_ind = ii;
           break
   elseif abs(y(ii,1)) < tol && start_flag == 1 && y(ii,1)<y(ii-1,1) && upward_flag == 0
            end_ind = ii;
            break
   end
end
T_y = end_ind - start_ind;
T_cycle = size(cycle,1);

if abs(T_y - T_cycle) > tol
    disp('new period');
end

phase_series = mod(0:(size(y,1)-start_ind),T_y)/T_y;

end

