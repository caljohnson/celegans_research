function [ updater_handle ] = discrete_neural_state_init_coupling(K_init, off_thresh, on_thresh, is_V)
%DISCRETE_NEURAL_STATE_INIT_COUPLING Initializes the discrete neural state updater
% returns a handle to a customized nested function 'discrete_neural_state_updater'.
% k_init specifies the initial value of the curvature, 
% thresh1, thresh2 define the thresholds of our neural state updater
% is_V =1 for V neurons, =0 for D neurons (swaps definition of handle)
%THIS ONE IS 50% SELF COUPLING, 50% NEIGHBOR COUPLING

%initialize state of neuron
if is_V == 1
    if( K_init >= on_thresh ) %if pass upper threshold, turn on
        S = 1;
    else %( K_init <= off_thresh ) %if pass lower threshold, turn off
        S = 0;
    end 
else %if is_V == 0, looking at D neuron, input is -K
    if( -K_init >= on_thresh ) %if pass upper threshold, turn on
        S = 1;
    else %( -K_init <= off_thresh ) %if pass lower threshold, turn off
        S = 0;
    end
end
    
%return updater handle, depending on if V or D neuron
if is_V == 1
    updater_handle = @(K1,K2) discrete_neural_state_updater(K1, K2);
else
    updater_handle = @(K1,K2) discrete_neural_state_updater(-K1, -K2);
end

%compute neural state and update
    function state = discrete_neural_state_updater(K1,K2)
        if( (50*K1+50*K2)/100 > on_thresh ) %if pass upper threshold, turn on
            state = 1;
            S = state;
        elseif( (50*K1+50*K2)/100 < off_thresh ) %if pass lower threshold, turn off
            state = 0;
            S = state;
        else
            state = S;
        end 
    end

end

