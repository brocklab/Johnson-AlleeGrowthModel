function [Nmodel] = singleexpmodeld(p, N0, t)

% This function is given a growth rate parameter p (1) and an initial cell 
% number N0 and returns the corresponding model cell numbers
num_reps = length(N0);
ind_rep = find(t == 0);
l = length(t)/num_reps;

% model
% N(t) = N_0 * exp(g*t)

for j = 1:num_reps
Nmod(:,j) = N0(j).* exp(-(p(1))*t(((j-1)*l)+1:j*l));
end
Nmodel = reshape(Nmod,[length(t),1]);

end