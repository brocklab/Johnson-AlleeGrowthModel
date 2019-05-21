%% Load in selected mudatavec
% Want to calculate BIC and AIC value from the fitting
% Need to make one big matrix of mu and vardata with corresponding time vecs in order
% Load in raw data, not down sampled
clear all; close all; clc

mudatavec = [];
vardatavec = [];
timevec = [];
mudatavecall = [];
vardatavecall = [];
timevecall = [];
N0vec = [];
numsamps = 0;

S = load('../out/BTsumfit.mat');
BTsum= S.BTsum;
N0list = [ 2 4 10];


for j = 1:length(N0list)
    for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
        mudata = BTsum(i).mu_t(2:end);
        ind = 1:1:length(mudata);
        mudatasmall= mudata(ind);
    mudatavec = vertcat(mudatavec, mudatasmall);
        vardata = BTsum(i).var_t(2:end);
        vardatasmall = vardata(ind);
    vardatavec = vertcat(vardatavec, vardatasmall);
        time = BTsum(i).timevec(2:end)';
        timesmall = time(ind);
    timevec = vertcat(timevec, timesmall);
  
   
    vardatavecall = vertcat(vardatavecall, BTsum(i).var_t(2:end));
    mudatavecall = vertcat(mudatavecall, BTsum(i).mu_t(2:end));
    timevecall = vertcat(timevecall, BTsum(i).timevec(2:end)');
    numsamps = BTsum(i).num + numsamps;

end
    end
end

%% confirm it's same
figure;
for j = 1:length(N0list)
plot(timevecall, mudatavecall, 'r*')
hold on
plot(timevec, mudatavec, 'b.')
end
xlabel('time (hours)')
ylabel('cell numer')
title('down sampled data for fitting')

%% Perform fitting: fit all data to a single b & d
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 1;
tsamp = 0:4:332; % set sample time (if data is down sampled change accordingly
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);
tsamp = [0 tsampnew];
%% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
theta = [bguess,dguess];
Ninit = N0list;
N= numsamps;
[pbestbd, BICbd, AICbd,negLLfitbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit b-d')
legend boxoff
title(['Mean  in data fit to simple birth-death model, BIC= ', num2str(BICbd)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit b-d')
legend boxoff
title(['Variance in data fit to simple birth-death model, BIC= ', num2str(BICbd)])
xlim([ 0 332])
%% Profile likelihood of b-d model 
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 2;

 % Check fxn

        fbest = negLLfitbd;
        factor = 0.05;
        numpoints = 10;
        params = pbestbd;
        [profilesbd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesbd(:,2,1)<threshold);
     
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIbd(1,:) = [NaN, NaN];
        else
             CIbd(1,:) = [profilesbd(ilo1(1),1,1), profilesbd(ilo1(end),1,1)];
        end
        ilo2 = find(profilesbd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIbd(2,:)= [NaN NaN];
        else
             CIbd(2,:) = [profilesbd(ilo2(1),1,2), profilesbd(ilo2(end),1,2)];
        end

    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter

%% Plot profiles to view CI

 figure;
         subplot(1, 2,1)
        plot(profilesbd(:,1,1), profilesbd(:,2,1))
        hold on
        plot(pbestbd(1),negLLfitbd,'r*','LineWidth',2)
        plot([profilesbd(1,1,1) profilesbd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiles param values')
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIbd(1,:)),']'])

        subplot(1,2,2)
        plot(profilesbd(:,1,2), profilesbd(:,2,2))
        hold on
        plot(pbestbd(2),negLLfitbd,'r*','LineWidth',2)
        plot([profilesbd(1,1,2) profilesbd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles param values')
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIbd(2,:)),']'])
    

%% Perform fitting: fit all data to strong Allee model on birth
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 2;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 2;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAb, BICstrAb, AICstrAb,negLLfitstrAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model fr strong Allee on birth
modelcode = 2;

figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to strong Allee (b), BIC= ', num2str(BICstrAb)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to strong Allee(b) model, BIC= ', num2str(BICstrAb)])
xlim([ 0 332])
%% Profile likelihood of strong Allee on birth model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAb;
        factor = 0.05;
        numpoints = 10;
        params = pbeststrAb;
        [profilesstrAb] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesstrAb(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(1,:) = [NaN, NaN];
        else
             CIstrAb(1,:) = [profilesstrAb(ilo1(1),1,1), profilesstrAb(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAb(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(2,:)= [NaN NaN];
        else
             CIstrAb(2,:) = [profilesstrAb(ilo2(1),1,2), profilesstrAb(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAb(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(3,:)= [NaN NaN];
        else
             CIstrAb(3,:) = [profilesstrAb(ilo3(1),1,3), profilesstrAb(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot profiles of strong Allee on birth to view CI

 figure;
         subplot(1, 3,1)
        plot(profilesstrAb(:,1,1), profilesstrAb(:,2,1))
        hold on
        plot(pbeststrAb(1),negLLfitstrAb,'r*','LineWidth',2)
        plot([profilesstrAb(1,1,1) profilesstrAb(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
        xlim([profilesstrAb(1,1,1) profilesstrAb(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIstrAb(1,:)),']'])

        subplot(1,3,2)
        plot(profilesstrAb(:,1,2), profilesstrAb(:,2,2))
        hold on
        plot(pbeststrAb(2),negLLfitstrAb,'r*','LineWidth',2)
        plot([profilesstrAb(1,1,2) profilesstrAb(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profilesstrAb(1,1,2) profilesstrAb(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIstrAb(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAb(:,1,3), profilesstrAb(:,2,3))
        hold on
        plot(pbeststrAb(3),negLLfitstrAb,'r*','LineWidth',2)
        plot([profilesstrAb(1,1,3) profilesstrAb(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAb(1,1,3) profilesstrAb(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIstrAb(3,:)),']'])
        


%% Perform fitting: fit all data to strong Allee model on death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 3;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAd, BICstrAd, AICstrAd,negLLfitstrAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to strong Allee (d), BIC= ', num2str(BICstrAd)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to strong Allee(d) model, BIC= ', num2str(BICstrAd)])
xlim([ 0 332])

%% Profile likelihood of strong Allee on death model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAd;
        factor = 0.05;
        numpoints = 10;
        params = pbeststrAd;
        [profilesstrAd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesstrAd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAd(1,:) = [NaN, NaN];
        else
             CIstrAd(1,:) = [profilesstrAd(ilo1(1),1,1), profilesstrAd(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAd(2,:)= [NaN NaN];
        else
             CIstrAd(2,:) = [profilesstrAd(ilo2(1),1,2), profilesstrAd(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAd(3,:)= [NaN NaN];
        else
             CIstrAd(3,:) = [profilesstrAd(ilo3(1),1,3), profilesstrAd(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter

%% Plot profiles of strong Allee on death to view CI

 figure;
         subplot(1, 3,1)
        plot(profilesstrAd(:,1,1), profilesstrAd(:,2,1))
        hold on
        plot(pbeststrAd(1),negLLfitstrAd,'r*','LineWidth',2)
        plot([profilesstrAd(1,1,1) profilesstrAd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIstrAd(1,:)),']'])

        subplot(1,3,2)
        plot(profilesstrAd(:,1,2), profilesstrAd(:,2,2))
        hold on
        plot(pbeststrAd(2),negLLfitstrAd,'r*','LineWidth',2)
        plot([profilesstrAd(1,1,2) profilesstrAd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profilesstrAd(1,1,2) profilesstrAd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIstrAd(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAd(:,1,3), profilesstrAd(:,2,3))
        hold on
        plot(pbeststrAd(3),negLLfitstrAd,'r*','LineWidth',2)
        plot([profilesstrAd(1,1,3) profilesstrAd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAd(1,1,3) profilesstrAd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIstrAd(3,:)),']'])
        



%% Perform fitting: fit all data to strong Allee model on  birth & death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 4;
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAbd, BICstrAbd, AICstrAbd,negLLfitstrAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to strong Allee ( b & d), BIC= ', num2str(BICstrAbd)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to strong Allee(b & d) model, BIC= ', num2str(BICstrAbd)])
xlim([ 0 332])


%% Profile likelihood of strong Allee on birth & death model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAbd;
        factor = 0.05;
        numpoints = 10;
        params = pbeststrAbd;
        [profilesstrAbd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesstrAbd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAbd(1,:) = [NaN, NaN];
        else
             CIstrAbd(1,:) = [profilesstrAbd(ilo1(1),1,1), profilesstrAbd(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAbd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAbd(2,:)= [NaN NaN];
        else
             CIstrAbd(2,:) = [profilesstrAbd(ilo2(1),1,2), profilesstrAbd(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAbd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAbd(3,:)= [NaN NaN];
        else
             CIstrAbd(3,:) = [profilesstrAbd(ilo3(1),1,3), profilesstrAbd(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot the profiles for strong Allee model on birth & death
    
    figure;
         subplot(1, 3,1)
        plot(profilesstrAbd(:,1,1), profilesstrAbd(:,2,1))
        hold on
        plot(pbeststrAbd(1),negLLfitstrAbd,'r*','LineWidth',2)
        plot([profilesstrAbd(1,1,1) profilesstrAbd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIstrAbd(1,:)),']'])

        subplot(1,3,2)
        plot(profilesstrAbd(:,1,2), profilesstrAbd(:,2,2))
        hold on
        plot(pbeststrAbd(2),negLLfitstrAbd,'r*','LineWidth',2)
        plot([profilesstrAbd(1,1,2) profilesstrAbd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        %xlim ([profilesstrAbd(1,1,2) profilesstrAbd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIstrAbd(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAbd(:,1,3), profilesstrAbd(:,2,3))
        hold on
        plot(pbeststrAbd(3),negLLfitstrAbd,'r*','LineWidth',2)
        plot([profilesstrAbd(1,1,3) profilesstrAbd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAbd(1,1,3) profilesstrAbd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIstrAbd(3,:)),']'])
        


%% Perform fitting: fit all data to weak/strong Allee on birth
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 5;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAb, BICwkAb, AICwkAb,negLLfitwkAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to weak Allee (b), BIC= ', num2str(BICwkAb)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to weak Allee(b) model, BIC= ', num2str(BICwkAb)])
xlim([ 0 332])

%% Profile likelihood of weak/strong Allee on birth 
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 4;

 % Check fxn

        fbest = negLLfitwkAb;
        factor = 0.05;
        numpoints = 10;
        params = pbestwkAb;
        [profileswkAb] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profileswkAb(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(1,:) = [NaN, NaN];
        else
             CIwkAb(1,:) = [profileswkAb(ilo1(1),1,1), profileswkAb(ilo1(end),1,1)];
        end
        ilo2 = find(profileswkAb(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(2,:)= [NaN NaN];
        else
             CIwkAb(2,:) = [profileswkAb(ilo2(1),1,2), profileswkAb(ilo2(end),1,2)];
        end
        ilo3 = find(profileswkAb(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(3,:)= [NaN NaN];
        else
             CIwkAb(3,:) = [profileswkAb(ilo3(1),1,3), profileswkAb(ilo3(end),1,3)];
        end
        ilo4 = find(profileswkAb(:,2,4)<threshold);
        if ilo4(end) ==numpoints*2 ||ilo4(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(4,:)= [NaN NaN];
        else
             CIwkAb(4,:) = [profileswkAb(ilo4(1),1,4), profileswkAb(ilo4(end),1,4)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot profiles for weak/strong Allee on birth
    
    figure;
         subplot(1, 4,1)
        plot(profileswkAb(:,1,1), profileswkAb(:,2,1))
        hold on
        plot(pbestwkAb(1),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,1) profileswkAb(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIwkAb(1,:)),']'])

        subplot(1,4,2)
        plot(profileswkAb(:,1,2), profileswkAb(:,2,2))
        hold on
        plot(pbestwkAb(2),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,2) profileswkAb(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profileswkAb(1,1,2) profileswkAb(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIwkAb(2,:)),']'])
        
        subplot(1,4,3)
        plot(profileswkAb(:,1,3), profileswkAb(:,2,3))
        hold on
        plot(pbestwkAb(3),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,3) profileswkAb(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profileswkAb(1,1,3) profileswkAb(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIwkAb(3,:)),']'])
        
        subplot(1,4,4)
        plot(profileswkAb(:,1,4), profileswkAb(:,2,4))
        hold on
        plot(pbestwkAb(4),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,4) profileswkAb(end,1,4)],[threshold threshold],'r--')
        xlabel(' profiled tau values')
        xlim([profileswkAb(1,1,4) profileswkAb(end,1,4)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of \tau CI = [',num2str(CIwkAb(4,:)),']'])



%% Perform fitting: fit all data to weak/strong Allee on death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 6;
tsamp = 0:4:332; % set sample time (if data is down sampled change accordingly
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);
tsamp = [0 tsampnew];
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAd, BICwkAd, AICwkAd,negLLfitwkAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to weak Allee (d), BIC= ', num2str(BICwkAd)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to weak Allee(d) model, BIC= ', num2str(BICwkAd)])
xlim([ 0 332])
%% Profile likelihood of weak/strong Allee on death
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 4;

 % Check fxn

        fbest = negLLfitwkAd;
        factor = 0.05;
        numpoints = 10;
        params = pbestwkAd;
        [profileswkAd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profileswkAd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(1,:) = [NaN, NaN];
        else
             CIwkAd(1,:) = [profileswkAd(ilo1(1),1,1), profileswkAd(ilo1(end),1,1)];
        end
        ilo2 = find(profileswkAd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(2,:)= [NaN NaN];
        else
             CIwkAd(2,:) = [profileswkAd(ilo2(1),1,2), profileswkAd(ilo2(end),1,2)];
        end
        ilo3 = find(profileswkAd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(3,:)= [NaN NaN];
        else
             CIwkAd(3,:) = [profileswkAd(ilo3(1),1,3), profileswkAd(ilo3(end),1,3)];
        end
        ilo4 = find(profileswkAd(:,2,4)<threshold);
        if ilo4(end) ==numpoints*2 ||ilo4(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(4,:)= [NaN NaN];
        else
             CIwkAd(4,:) = [profileswkAd(ilo4(1),1,4), profileswkAd(ilo4(end),1,4)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
    %%
    save('../out/profilesbd.mat', 'profilesbd')
    save('../out/profilesstrAb.mat', 'profilesstrAb')
    save('../out/profilesstrAd.mat', 'profilesstrAd')
    save('../out/profilesstrAbd.mat', 'profilesstrAbd')
    save('../out/profileswkAb.mat', 'profileswkAb')
    save('../out/profileswkAd.mat', 'profileswkAd')
%% Plot profiles for weak/strong Allee on death
    
    figure;
         subplot(1, 4,1)
        plot(profileswkAd(:,1,1), profileswkAd(:,2,1))
        hold on
        plot(pbestwkAd(1),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,1) profileswkAd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIwkAd(1,:)),']'])

        subplot(1,4,2)
        plot(profileswkAd(:,1,2), profileswkAd(:,2,2))
        hold on
        plot(pbestwkAd(2),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,2) profileswkAd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        %xlim ([profileswkAd(1,1,2) profileswkAd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIwkAd(2,:)),']'])
        
        subplot(1,4,3)
        plot(profileswkAd(:,1,3), profileswkAd(:,2,3))
        hold on
        plot(pbestwkAd(3),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,3) profileswkAd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profileswkAd(1,1,3) profileswkAd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIwkAd(3,:)),']'])
        
        subplot(1,4,4)
        plot(profileswkAd(:,1,4), profileswkAd(:,2,4))
        hold on
        plot(pbestwkAd(4),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,4) profileswkAd(end,1,4)],[threshold threshold],'r--')
        xlabel(' profiled tau values')
        xlim([profileswkAd(1,1,4) profileswkAd(end,1,4)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of \tau CI = [',num2str(CIwkAd(4,:)),']'])


%% Perform fitting: fit all data to weak/strong Allee on birth & death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 7;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAbd, BICwkAbd, AICwkAbd,negLLfitwkAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to weak Allee (b & d), BIC= ', num2str(BICwkAbd)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to weak Allee(b & d) model, BIC= ', num2str(BICwkAbd)])
xlim([ 0 332])

%% Profile likelihood of weak/strong Allee on birth & death
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 4;

 % Check fxn

        fbest = negLLfitwkAbd;
        factor = 0.05;
        numpoints = 10;
        params = pbestwkAbd;
        [profileswkAbd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profileswkAbd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(1,:) = [NaN, NaN];
        else
             CIwkAbd(1,:) = [profileswkAbd(ilo1(1),1,1), profileswkAbd(ilo1(end),1,1)];
        end
        ilo2 = find(profileswkAbd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(2,:)= [NaN NaN];
        else
             CIwkAbd(2,:) = [profileswkAbd(ilo2(1),1,2), profileswkAbd(ilo2(end),1,2)];
        end
        ilo3 = find(profileswkAbd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(3,:)= [NaN NaN];
        else
             CIwkAbd(3,:) = [profileswkAbd(ilo3(1),1,3), profileswkAd(ilo3(end),1,3)];
        end
        ilo4 = find(profileswkAbd(:,2,4)<threshold);
        if ilo4(end) ==numpoints*2 ||ilo4(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(4,:)= [NaN NaN];
        else
             CIwkAbd(4,:) = [profileswkAbd(ilo4(1),1,4), profileswkAbd(ilo4(end),1,4)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter


%% Plot profiles for weak/strong Allee on birth death
    
    figure;
         subplot(1, 4,1)
        plot(profileswkAbd(:,1,1), profileswkAbd(:,2,1))
        hold on
        plot(pbestwkAbd(1),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,1) profileswkAbd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIwkAbd(1,:)),']'])

        subplot(1,4,2)
        plot(profileswkAbd(:,1,2), profileswkAbd(:,2,2))
        hold on
        plot(pbestwkAbd(2),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,2) profileswkAbd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        %xlim ([profileswkAd(1,1,2) profileswkAd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIwkAbd(2,:)),']'])
        
        subplot(1,4,3)
        plot(profileswkAbd(:,1,3), profileswkAbd(:,2,3))
        hold on
        plot(pbestwkAbd(3),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,3) profileswkAbd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profileswkAbd(1,1,3) profileswkAbd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIwkAbd(3,:)),']'])
        
        subplot(1,4,4)
        plot(profileswkAbd(:,1,4), profileswkAbd(:,2,4))
        hold on
        plot(pbestwkAbd(4),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,4) profileswkAbd(end,1,4)],[threshold threshold],'r--')
        xlabel(' profiled tau values')
        xlim([profileswkAbd(1,1,4) profileswkAbd(end,1,4)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of \tau CI = [',num2str(CIwkAbd(4,:)),']'])

   save('../out/profileswkAbd.mat', 'profileswkAbd')

%% Make vector of BIC values, negLL values, and parameter numbers
num_params = [ 2, 3, 3, 3, 4, 4, 4];
modellist = [1, 2, 3, 4, 5, 6, 7];
BICvals = [ BICbd, BICstrAb, BICstrAd, BICstrAbd, BICwkAb, BICwkAd, BICwkAbd];
negLLvals = [negLLfitbd, negLLfitstrAb, negLLfitstrAd, negLLfitstrAbd, negLLfitwkAb, negLLfitwkAd, negLLfitwkAbd];
modelnames = {'  ' 'b-d', 'strAb', 'strAd', 'strAbd', 'wkAb', 'wkAd', 'wkAbd', '   '};
shortnames = { 'strAb', 'wkAb'};
figure;
subplot(2,1,1)
semilogy(modellist, BICvals, 'ro', 'LineWidth', 4)
xlabel('Model')
ylabel('BIC value')
xlim([0 8])
ylim([ min(BICvals)-300 max(BICvals)+300])
title('BIC values for each model (sampling every 16 hours)')
set(gca,'xticklabel',modelnames)


% subplot(3,1,2)
% semilogy(modellist, negLLvals, 'ro', 'LineWidth', 4)
% xlabel('Model')
% ylabel(' best negLL')
% xlim([0 8])
% ylim([ min(negLLvals)-5, max(negLLvals)+5])
% title('negLL values for each model (sampling every 4 hours)')
% set(gca,'xticklabel',modelnames)

subplot(2,1,2)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
ylabel(' number of parameters')
xlim([0 8])
ylim([ 1 5 ])
title ('number of parameters')
set(gca,'xticklabel',modelnames)
figure;
plot([1 2], [(BICvals(2)) (BICvals(5))], 'ro', 'LineWidth',4)
xlabel('Model')
% xlim([ 0 3])
ylim([ 1800 1860])
ylabel('BIC value')
title('BIC model comparison')
set(gca,'xticklabel',shortnames)

