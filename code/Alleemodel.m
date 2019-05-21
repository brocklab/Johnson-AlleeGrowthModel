function [ Nmodellong] = Alleemodel(P, tlong, N0)
%Alleemodel This function outputs a single vector that corresponds to the
%long time vector given to it

%   Detailed explanation goes here
    % inputs = params of Allee model, long time vector, initial cell
    % numbers
    P = num2cell(P);
[g, A, carcap] = deal(P{:}); % our parameters

num_reps = length(N0);

Nmodellong = [];
t0 = find(ismember(tlong, 0)); % find start of new sample in time vector

for i = 1:num_reps-1
  t = tlong(t0(i):(t0(i+1)-1)); % only grab from new sample and input corresponding N0
  [tout,Nmodel] = ode45(@(t,Nmodel) odefunAllee(t, Nmodel,g, A, carcap ),t, N0(i));
   Nmodellong = vertcat(Nmodellong, Nmodel);
end
for i = num_reps
   t = tlong(t0(i):end);
   [tout,Nmodel] = ode45(@(t,Nmodel) odefunAllee(t, Nmodel,g, A, carcap ),t, N0(i));
   Nmodellong = vertcat(Nmodellong, Nmodel);
end

% Remove any negative values ( model can't go negative)
for j = 1:length(Nmodellong)
    if Nmodellong(j) <0
        Nmodellong(j) = 0;
    end
end

end




