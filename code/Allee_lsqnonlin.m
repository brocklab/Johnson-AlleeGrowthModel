%% Allee effect fitting
% Use Nfake and see if lsqnonlin can capture input parameters of N0, g, and
% A
LB = zeros(2,1);  % Lower Bounds
UB = [Inf Inf]; % Upper Bounds
params0 = [.01; 7];% Initial Guess: N0, g, A 
% Actual params: g = 0.01 A = 7
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);


paramsAllee= lsqnonlin(@fit_Allee_closed, params0, LB, UB, options, Nfake, N0, tbig);
Afit = paramsAllee(2);
gfit = paramsAllee(1);
for i = 1:length(N0)
Nfit(:,i)= (N0(i)-Afit).*exp(gfit*tbig(i,:)) + Afit;
for j = 1:length(Nfit(:,i))
    if Nfit(j,i) <0
        Nfit(j,i) = 0;
    end
end
end

figure;
hold off
for j = 1:length(N0)
    hold on
    semilogy(tbig(j,:), Nfit(:,j), 'LineWidth', 3)
    hold on
    semilogy(tbig(j,:), Nfake(:,j), 'LineWidth', 2)
    text(tbig(j, end-27), Nfit(end-27,j), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
   % plot(tout, Nfake(:,j), 'o')
   %ylim([ 0 200])
end
xlabel('time')
ylabel('N')
title(['Fit growth dynamics versus input data, A_{fit}= ', num2str(Afit), ', A_{true}= ', num2str(A), ', \eta=', num2str(eta)])
%legend('N_{0} = 2', 'N_{0} = 5', 'N_{0} = 10', 'N_{0} = 12', 'N_{0} = 20')
legend boxoff
