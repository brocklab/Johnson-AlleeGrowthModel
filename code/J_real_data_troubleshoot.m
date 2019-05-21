% Trouble shooting J function with fitting for sigma
% this code is meant to go through each piece of the negative LL function
% and trouble shoot it, using the neg_LL as a function of time for the
% guessed parameters, and incorrect but converged parameters.

% Run simultaneously with moment_approach_fitting.m

% current J function
J4_t= @(phat)((0.5*((log(2*pi*var_in_mean4(pbxform4(phat)))) + (((yfxform(modelfun_mu4(pbxform4(phat))) - yfxform(mudata))./sqrt(var_in_mean4(pbxform4(phat)))).^2))) +...
 (0.5*((log(2*pi*var_in_var4(pbxform4(phat)))) + (((yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata))./sqrt(var_in_var4(pbxform4(phat)))).^2))));

%%
guess = [];
fit = [];
% Look only at mean, model - data
J1_4 = @(phat)(((yfxform(modelfun_mu4(pbxform4(phat))) - yfxform(mudata)).^2)./(var_in_mean4(pbxform4(phat))));
guess(:,1) = J1_4(pfxform4(theta4));
fit(:,1) = J1_4(phatbest_J4);

figure;
plot(t, mudata, 'LineWidth',2)
hold on
plot(t, modelfun_mu4(pbxform4(phatbest_J4)), 'LineWidth',2)
legend('data', 'fit params')


figure;
hold on
plot(t,guess(:,1),'*','LineWidth', 2)
plot(t,fit(:,1), 'LineWidth', 2)
legend( 'initial guess', 'fit params')
xlabel('time')
ylabel('component 1 of neg LL')
title('squared error/var in mean')
%%
J2_4 = @(phat)((log(2*pi*(var_in_mean4(pbxform4(phat)))))+ ((yfxform(modelfun_mu4(pbxform4(phat))) - yfxform(mudata)).^2)./(var_in_mean4(pbxform4(phat))));
guess(:,2) = J2_4(pfxform4(theta4));
fit(:,2) = J2_4(phatbest_J4);

figure;
hold on
plot(t,guess(:,2),'*','LineWidth', 2)
plot(t,fit(:,2), 'LineWidth', 2)
legend( 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in mean + squared error/var in mean')

J3_4 = @(phat)((log(2*pi*(var_in_mean4(pbxform4(phat))))));
guess(:,3) = J3_4(pfxform4(theta4));
fit(:,3) = J3_4(phatbest_J4);
new= J3_4(pfxform4(horzcat(params_best_J4(1:2),theta4(3:4))));

figure;
hold on
plot(t,guess(:,3),'*','LineWidth', 2)
plot(t,fit(:,3), 'LineWidth', 2)
plot(t, new, 'LineWidth',2)
legend( 'initial guess', 'fit params', 'test')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in mean')

%% Next focus on variance terms

J4_4 = @(phat) ((yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata)).^2)./(var_in_var4(pbxform4(phat)));

guess(:,4) = J4_4(pfxform4(theta4));
fit(:,4) = J4_4(phatbest_J4);

figure;
hold on
plot(t,guess(:,4),'*','LineWidth', 2)
plot(t,fit(:,4), 'LineWidth', 2)
legend('initial guess', 'fit params')
xlabel('time')
ylabel('component of neg LL')
xlim( [0 100])
title('squared error in variance/var in var')


%%
J5_4 = @(phat) (log(2*pi*(var_in_var4(pbxform4(phat))))+(yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata)).^2./(var_in_var4(pbxform4(phat))));
guess(:,5) = J5_4(pfxform4(theta4));
fit(:,5) = J5_4(phatbest_J4);

figure;
hold on
plot(t,guess(:,5),'*','LineWidth', 2)
plot(t,fit(:,5), 'LineWidth', 2)
legend('initial guess', 'fit params')
xlabel('time')
ylabel('component of neg LL')
title(' log 2pi var in var + squared error in variance/var in var')

J7_4 = @(phat)log(2*pi*(var_in_var4(pbxform4(phat))))
guess(:,6) = J7_4(pfxform4(theta4));
fit(:,6) = J7_4(phatbest_J4);

figure;
hold on
plot(t,guess(:,6),'*','LineWidth', 2)
plot(t,fit(:,6), 'LineWidth', 2)
legend('initial guess', 'fit params')
xlabel('time')
ylabel('component of neg LL')
title(' log 2pi var in var')
%%
J6_4 = @(phat) ((yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata)).^2);
guess(:,6) = J6_4(pfxform4(theta4));
fit(:,6) = J6_4(phatbest_J4);

figure;
hold on
plot(t,guess(:,6),'*','LineWidth', 2)
plot(t,fit(:,6), 'LineWidth', 2)
legend('initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('squared error in variance')
%%
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

