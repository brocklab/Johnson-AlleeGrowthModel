% Trouble shooting J function
% this code is meant to go through each piece of the negative LL function
% and trouble shoot it, using the neg_LL as a function of time for the true
% parameters, guessed parameters, and incorrect but converged parameters.

% Run simultaneously with moment_approach_fitting.m

% current J function
J_t= @(phat)((0.5*((log(2*pi*var_in_mean(pbxform(phat)))) + (((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
 (0.5*((log(2*pi*var_in_var(pbxform(phat)))) + (((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata))./sqrt(var_in_var(pbxform(phat)))).^2))));
% current J_t (no sum)
J_t= @(phat)((0.5*((log(2*pi*var_in_mean(pbxform(phat)))) + (((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
 (0.5*((log(2*pi*var_in_var(pbxform(phat)))) + (((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata))./sqrt(var_in_var(pbxform(phat)))).^2))));
%%
guess = [];
true = [];
fit = [];
% Look only at mean, model - data
J1 = @(phat)(((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata)).^2)./(var_in_mean(pbxform(phat))));
true(:,1) = J1(pfxform(p));
guess(:,1) = J1(pfxform(theta));
fit(:,1) = J1(phatbest_J);

figure;
plot(t, mudata)
hold on
plot(t,modelfun_mu(p))
plot(t, modelfun_mu(theta))


figure;
plot(t, true(:,1), 'LineWidth',2)
hold on
plot(t,guess(:,1),'*','LineWidth', 2)
plot(t,fit(:,1), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('component 1 of neg LL')
title('squared error/var in mean')

J2 = @(phat)((log(2*pi*(var_in_mean(pbxform(phat)))))+ ((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata)).^2)./(var_in_mean(pbxform(phat))));
true(:,2) = J2(pfxform(p));
guess(:,2) = J2(pfxform(theta));
fit(:,2) = J2(phatbest_J);

figure;
plot(t, true(:,2), 'LineWidth',2)
hold on
plot(t,guess(:,2),'*','LineWidth', 2)
plot(t,fit(:,2), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in mean + squared error/var in mean')

J3 = @(phat)((log(2*pi*(var_in_mean(pbxform(phat))))));
true(:,3) = J3(pfxform(p));
guess(:,3) = J3(pfxform(theta));
fit(:,3) = J3(phatbest_J);

figure;
plot(t, true(:,3), 'LineWidth',2)
hold on
plot(t,guess(:,3),'*','LineWidth', 2)
plot(t,fit(:,3), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in mean')

%% Next focus on variance terms

J4 = @(phat) ((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata)).^2)./(var_in_var(pbxform(phat)));
true(:,4) = J4(pfxform(p));
guess(:,4) = J4(pfxform(theta));
fit(:,4) = J4(phatbest_J);

figure;
plot(t, true(:,4), 'LineWidth',2)
hold on
plot(t,guess(:,4),'*','LineWidth', 2)
plot(t,fit(:,4), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('component of neg LL')
xlim( [0 100])
title('squared error in variance/var in var')


%%
J5 = @(phat) (log(2*pi*(var_in_var(pbxform(phat))))+(yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata)).^2./(var_in_var(pbxform(phat))));
true(:,5) = J5(pfxform(p));
guess(:,5) = J5(pfxform(theta));
fit(:,5) = J5(phatbest_J);

figure;
plot(t, true(:,5), 'LineWidth',2)
hold on
plot(t,guess(:,5),'*','LineWidth', 2)
plot(t,fit(:,5), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('component of neg LL')
title(' log 2pi var in var + squared error in variance/var in var')
%%
J6 = @(phat)((log(2*pi*(var_in_var(pbxform(phat))))));
true(:,6) = J6(pfxform(p));
guess(:,6) = J6(pfxform(theta));
fit(:,6) = J6(phatbest_J);

figure;
plot(t, true(:,6), 'LineWidth',2)
hold on
plot(t,guess(:,6),'*','LineWidth', 2)
plot(t,fit(:,6), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in var')
%%
true(:,7) = var_in_var(p);
guess(:,7) = var_in_var(theta);
fit(:,7) = var_in_var(params_best_J);

figure;
plot(t, true(:,7), 'LineWidth', 2)
hold on
plot(t, var_in_var_data, 'LineWidth',2)
legend('var in var fxn', 'var in var data')
xlabel('time')
ylabel('var in var')
title('Comparison data to function')
%%

figure;
plot(t, true(:,7), 'LineWidth',2)
hold on
plot(t,guess(:,7),'*','LineWidth', 2)
%plot(t,fit(:,7), 'LineWidth', 2)
plot(t, var_in_var_data, 'o','LineWidth', 2)
legend('true params', 'initial guess', 'var in var data')
xlabel('time')
ylabel('var in var')
title('var in var')
%%
figure;
plot(t, 1/true(:,7), 'LineWidth',2)
hold on
plot(t,1/guess(:,7),'*','LineWidth', 2)
plot(t,1/fit(:,7), 'LineWidth', 2)
legend('true params', 'initial guess', 'fit params')
xlabel('time')
ylabel('1/var in var')
title('1/var in var')

