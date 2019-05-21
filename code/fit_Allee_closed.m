function [ err_long ] = fit_Allee_closed( params, Nmeas, N0, tbig)
%This function returns a cost function to me minimized
% Nmeas is a matrix that contains number of sample rows by the number of
% time point columns
g = params(1);
A = params(2);

% Need run through a loop for each sample to produce curve

for i = 1:length(N0)
    % isolate each sample to be run through model
    t = tbig(i,:);
    Nmodel(:,i) = (N0(i)-A).*exp(g*t) + A;
    
end



% now reshape both of them into vector form

 
 err = Nmodel - Nmeas;
 err_long = reshape(err, (size(err,1)*size(err,2)),1);
 



end