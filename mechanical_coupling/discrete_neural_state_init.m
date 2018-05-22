function [ updater_handle ] = discrete_neural_state_init( S_init, K_init, off_thresh, on_thresh, is_V)
%DISCRETE_NEURAL_STATE_INIT Initializes the discrete neural state updater
% returns a handle to a customized nested function 'discrete_neural_state_updater'.
% k_init specifies the initial value of the curvature, 
% thresh1, thresh2 define the thresholds of our neural state updater
% is_V =1 for V neurons, =0 for D neurons (swaps definition of handle)

%initialize state of neuron
if( K_init > on_thresh ) %if pass upper threshold, turn on
    S = 1;
elseif( K_init < off_thresh ) %if pass lower threshold, turn off
    S = 0;
else
    S = S_init;
end 

%return updater handle, depending on if V or D neuron
if is_V == 1
    updater_handle = @(K) discrete_neural_state_updater(K);
else
    updater_handle = @(K) discrete_neural_state_updater(-K);
end

%compute neural state and update
    function state = discrete_neural_state_updater(K)
        if( K > on_thresh ) %if pass upper threshold, turn on
            state = 1;
            S = state;
        elseif( K < off_thresh ) %if pass lower threshold, turn off
            state = 0;
            S = state;
        else
            state = S;
        end 
    end

end

