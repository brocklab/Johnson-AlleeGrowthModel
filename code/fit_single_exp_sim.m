function [diff] = fit_single_exp_sim(g, state, tstate)


% Here state contains 100 columns of discrete cell trajectories
% tstate contains 100 columns of times corresponding to those trajectories
num_sims = size(tstate,2);
num_iters = size(tstate,1);
vert_length = num_sims*num_iters;
for j =1:num_sims
    N0 = state(1,j);
    tmeas = tstate(:,j);
    N_model(:,j) = singleexpmodel(g, state(1,j),tmeas);
end

Nmeas = reshape(state, vert_length,1);
N_modellong = reshape(N_model, vert_length,1);
diff = N_modellong - Nmeas;

end

