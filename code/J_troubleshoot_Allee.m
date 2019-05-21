% Trouble shooting J function for Allee model
% this code is meant to go through each piece of the negative LL function
% and trouble shoot it, using the neg_LL as a function of time for the true
% parameters, guessed parameters, and incorrect but converged parameters.

% Run simultaneously with moment_approach_fitting.m

% current J function
J_t= @(phat)((0.5*((log(2*pi*var_in_mean5(pbxform5(phat)))) + (((yfxform(modelfun_mu5(pbxform5(phat))) - yfxform(mudata))./sqrt(var_in_mean5(pbxform5(phat)))).^2))) +...
 (0.5*((log(2*pi*var_in_var5(pbxform5(phat)))) + (((yfxform(modelfun_V5(pbxform5(phat))) - yfxform(vardata))./sqrt(var_in_var5(pbxform5(phat)))).^2))));
%%
guess = [];
fit = [];
% Look only at mean, model - data
%J1_5= @(phat)(((yfxform(modelfun_mu5(pbxform5(phat))) - yfxform(mudata))./sqrt(var_in_mean5(pbxform5(phat)))).^2); 
J1_5 = @(phat)((((modelfun_mu5(pbxform5(phat))) - (mudata)).^2)./(var_in_mean5(pbxform5(phat))));
%J1_5 = @(phat)(((modelfun_mu5(pbxform5(phat))) - (mudata)).^2);
guess(:,1) = J1_5(pfxform5(theta5));
fit(:,1) = J1_5(phatbest_Jnew);
true (:,1) = J1_5(pfxform5(p5));


g = sum(guess(:,1))
f = sum(fit(:,1))
tr = sum(true(:,1))
figure;
hold on
%%plot(t,guess(:,1),'*','LineWidth', 2)
plot(t,fit(:,1), 'LineWidth', 2)
plot(t, true (:,1), '*','LineWidth',2)
legend( 'fit params', 'true params')
xlabel('time')
ylabel('component 1 of neg LL')
title('squared error/var in mean')

figure;
hold on
plot(t, mudata, '*')
plot(t, modelfun_mu5(p5), '-')
plot(t, modelfun_mu5(pbxform5(phatbest_Jnew)), '-')
legend ('data', 'model true params', 'model fit params')
%%
J2_5= @(phat)(0.5.*((log(2*pi*var_in_mean5(pbxform5(phat)))) + ((((modelfun_mu5(pbxform5(phat))) - (mudata)).^2)./(var_in_mean5(pbxform5(phat))))));
guess(:,2) = J2_5(pfxform5(theta5));
fit(:,2) = J2_5(phatbest_Jnew);
true(:,2) = J2_5(pfxform5(p5));

figure;
hold on
plot(t,guess(:,2),'*','LineWidth', 2)
plot(t,fit(:,2), 'LineWidth', 2)
plot(t, true(:,2),'*', 'LineWidth',2)
legend( 'initial guess', 'fit params', 'true params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in mean + squared error/var in mean')
%%
J3_5 = @(phat)((log(2*pi*(var_in_mean5(pbxform5(phat))))));
guess(:,3) = J3_5(pfxform5(theta5));
fit(:,3) = J3_5(phatbest_J5);

figure;
hold on
plot(t,guess(:,3),'*','LineWidth', 2)
plot(t,fit(:,3), 'LineWidth', 2)
legend( 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in mean')
%% Check variance
J4_5 = @(phat)(((yfxform(modelfun_V5(pbxform5(phat))) - yfxform(vardata)).^2)./(var_in_var5(pbxform5(phat))));
guess(:,4) = J4_5(pfxform5(theta5));
fit(:,4) = J4_5(phatbest_J5);

figure;
hold on
plot(t,guess(:,4),'*','LineWidth', 2)
plot(t,fit(:,4), 'LineWidth', 2)
legend( 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('squared error/var in var')

J5_5= @(phat)((0.5*((log(2*pi*var_in_var5(pbxform5(phat)))) + (((yfxform(modelfun_V5(pbxform5(phat))) - yfxform(vardata))./sqrt(var_in_var5(pbxform5(phat)))).^2))));
guess(:,5) = J5_5(pfxform5(theta5));
fit(:,5) = J5_5(phatbest_J5);

figure;
hold on
plot(t,guess(:,5),'*','LineWidth', 2)
plot(t,fit(:,5), 'LineWidth', 2)
legend( 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in var + squared error/var in var')

J6_5= @(phat)((0.5*((log(2*pi*var_in_var5(pbxform5(phat)))))));
guess(:,6) = J6_5(pfxform5(theta5));
fit(:,6) = J6_5(phatbest_J5);

figure;
hold on
plot(t,guess(:,6),'*','LineWidth', 2)
plot(t,fit(:,6), 'LineWidth', 2)
legend( 'initial guess', 'fit params')
xlabel('time')
ylabel('component 2 of neg LL')
title('log 2pi*var in var')