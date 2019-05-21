% Process individual well data and wells grouped by Nseed data

% load fitted and processed data
close all; clear all; clc
S = load('../out/greenBT.mat');
green= S.green;
S = load('../out/greensumBT.mat');
greensum= S.greensum;

 %%       
figure (1)
subplot(2,6,1:6)
for j = 1:length(green)
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('All wells')

subplot(2,6,7)
for j = 1:length(green)
    if green(j).Nseed ==4
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
ylim([0, 1000])
xlabel('Time (hours)')
ylabel('N')
title('N_{seed} = 4')

subplot(2,6,8)
for j = 1:length(green)
    if green(j).Nseed ==8
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('N_{seed} = 8')


subplot(2,6,9)
for j = 1:length(green)
    if green(j).Nseed ==16
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('N_{seed} = 16')

subplot(2,6,10)
for j = 1:length(green)
    if green(j).Nseed ==32
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
ylim([0, 1000])
xlabel('Time (hours)')
ylabel('N')
title('N_{seed} = 32')

subplot(2,6,11)
for j = 1:length(green)
    if green(j).Nseed ==64
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('N_{seed} = 64')


subplot(2,6,12)
for j = 1:length(green)
    if green(j).Nseed ==128
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('N_{seed} = 128')
%%
figure;
hold off
subplot(1,3,1)
for j = 1:length(green)
    time = green(j).time;
    if green(j).takeoff == 1
    plot(green(j).time, log(green(j).cellnum/green(j).Nseed), 'color', green(j).cc, 'LineWidth', 1)
    end
    hold on
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Raw data takeoffs only')

subplot(1,3,2)
for j = 1:length(green)
    time = green(j).time;
    if green(j).takeoff == 1
    plot(green(j).time, log(green(j).Nmodellsq/green(j).Nseed), 'color', green(j).cc, 'LineWidth', 1)
    end
    hold on
end
min_time =green(j).time(end);
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Per capita growth takeoffs')



subplot(1,3,3)
hold off
for j = 1:length(greensum)
    model_low =  singleexpmodel(greensum(j).mean_gtkoff - 1.96*sqrt(greensum(j).var_gtkoff),greensum(j).Nseed, time);
    model_high =  singleexpmodel(greensum(j).mean_gtkoff + 1.96*sqrt(greensum(j).var_gtkoff),greensum(j).Nseed, time);
    plot(time, greensum(j).model_gtkoff, 'color', greensum(j).cc, 'LineWidth', 3)

    hold on
    %plot(time, model_low, 'color', greensum(j).cc, 'LineWidth', 1)
    %plot(time, model_high, 'color', greensum(j).cc, 'LineWidth', 1)
   %text(time(end), log(greensum(j).model_gtkoff(end)./greensum(j).Nseed), ['N_{0}= ', num2str(greensum(j).Nseed), ', g_{mean}= ', num2str(round(greensum(j).mean_gtkoff,4))], 'HorizontalAlignment', 'Right')
    hold on
end
legend('N_{seed}= 8', 'N_{seed}= 16', 'N_{seed}= 32')
legend boxoff
xlim([0, min_time])
xlabel('Time (hours)')
ylabel('N')
title('Mean population growth dynamics')
%Put error bars on growth rate estimates
%%
tint = [42 80 120];
for i = 1:3
ind(i) = find(ismember( round(time), tint(i)));
end
%%
subplot(1,3,3)
hold off
for j = 1:length(greensum)
    model_low =  singleexpmodel(greensum(j).mean_gtkoff - sqrt(greensum(j).var_gtkoff),greensum(j).Nseed, time);
    model_high =  singleexpmodel(greensum(j).mean_gtkoff + sqrt(greensum(j).var_gtkoff),greensum(j).Nseed, time);
    plot(time, log(greensum(j).model_gtkoff./greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 3)
    hold on
end
for j = 1:length(greensum)
    errorbar([42;80;120], log(greensum(j).model_gtkoff(ind)./greensum(j).Nseed), 0.5.*(log(model_high(ind)./greensum(j).Nseed)-log(model_low(ind)./greensum(j).Nseed)), 'color', greensum(j).cc, 'LineWidth', 1.3)
    %plot(time, log(model_low/greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 1)
    %plot(time, log(model_high/greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 1)
    %text(time(end), log(greensum(j).model_gtkoff(end)./greensum(j).Nseed), ['N_{0}= ', num2str(greensum(j).Nseed), ', g_{mean}= ', num2str(round(greensum(j).mean_gtkoff,4)),' +/- ', num2str(round(sqrt(greensum(j).var_gtkoff),3)),' cells/hr'], 'HorizontalAlignment', 'Right')
    hold on
end
legend('N_{seed}= 8', 'N_{seed}= 16', 'N_{seed}= 32')
legend boxoff
xlim([0, min_time])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Average Per Capita Growth by N_{seed}')

for j = 1:length(greensum)
    norm_var(j) = greensum(j).norm_var;
end

figure;
subplot(3,3,[2,5,9])
hold off
for j = 1:length(greensum)
    model_low =  singleexpmodel(greensum(j).mean_gtkoff - sqrt(greensum(j).var_gtkoff),greensum(j).Nseed, time);
    model_high =  singleexpmodel(greensum(j).mean_gtkoff + sqrt(greensum(j).var_gtkoff),greensum(j).Nseed, time);
    plot(time, log(greensum(j).model_gtkoff./greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 3)
    hold on
end
for j = 1:length(greensum)
    errorbar([42;80;120], log(greensum(j).model_gtkoff(ind)./greensum(j).Nseed), 0.5.*(log(model_high(ind)./greensum(j).Nseed)-log(model_low(ind)./greensum(j).Nseed)), 'color', greensum(j).cc, 'LineWidth', 1.3)
    %plot(time, log(model_low/greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 1)
    %plot(time, log(model_high/greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 1)
    %text(time(end), log(greensum(j).model_gtkoff(end)./greensum(j).Nseed), ['N_{0}= ', num2str(greensum(j).Nseed), ', g_{mean}= ', num2str(round(greensum(j).mean_gtkoff,4)),' +/- ', num2str(round(sqrt(greensum(j).var_gtkoff),3)),' cells/hr'], 'HorizontalAlignment', 'Right')
    hold on
end
legend('N_{seed}= 8', 'N_{seed}= 16', 'N_{seed}= 32')
legend boxoff
xlim([0, min_time])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Average Per Capita Growth by N_{seed}')
subplot(3,3,1)
for j = 1:length(green)
    if green(j).Nseed ==8
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
ylim([0, 1000])
xlabel('Time (hours)')
ylabel('N')
%title(['N_{seed} = 8'])%, Pct_{take-off}= ', num2str(round(greensum(1).pct_takeoff,0)), '%'])

subplot(3,3,4)
for j = 1:length(green)
    if green(j).Nseed ==16
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end
xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
%title(['N_{seed} = 16'])% Pct_{take-off}= ', num2str(round(greensum(2).pct_takeoff,0)), '%'])


subplot(3,3,7)
for j = 1:length(green)
    if green(j).Nseed ==32
    plot(green(j).time, green(j).cellnum, 'color', green(j).cc, 'LineWidth', 1)
    hold on
    end
end

xlim([0, green(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
%title(['N_{seed} = 32'])%, Pct_{take-off}= ', num2str(round(greensum(3).pct_takeoff, 0)), '%'])
%%
figure;
bar(norm_var)
ylabel('Normalized Variance')
%% Analyze 8 cell wells by actual initial seeding

num_8s = 30;
for j = 1:num_8s
    gs(j) = green(j).gtot;
    nums(j) = green(j).N0actual;
    cols(j) = green(j).num_colonies;
end

b1c = mldivide(nums', cols');
b1g = mldivide(nums',gs');
yfitc = b1c*nums;
yfitg = b1g*nums;
yresidc = cols - yfitc;
yresidg = gs-yfitg;
%Square the residuals and total them to obtain the residual sum of squares:
SSresidc = sum(yresidc.^2);
SSresidg = sum(yresidg.^2);
%Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotalc = (length(cols)-1) * var(cols);
SStotalg = (length(gs)-1)*var(gs);
%Compute R2 using the formula given in the introduction of this topic:
rsqc = 1 - (SSresidc./SStotalc);
rsqg = 1-(SSresidg./SStotalg);



subplot(1,2,1)
hold off
plot(nums, gs, 'g*')
hold on
plot(nums, yfitg, 'g-')
xlabel ('Actual N_{seed}')
ylabel('growth rate (cells/hr) of well')
title(['Does growth rate scale on N_{seed}? Rsq = ',num2str(rsqg)])
xlim ([0 8])

subplot(1,2,2)
hold off
plot(nums, cols', 'b*')
hold on
plot(nums, yfitc, '-')
xlabel ('Actual N_{seed}')
ylabel('Final number of colonies formed')
title(['Does number of colonies formed scale on N_{seed}? Rsq = ',num2str(rsqc)])
xlim ([0 8])

[val, count] = count_unique(nums);
mu = mean(nums);
sigma = std(nums);

figure;
hold off
bar(val, count)
xlabel('Actual number of cells seeded')
ylabel('Frequency')
title(['Average and standard deviation of actual cells seeded when N_{0}=8, \mu=,' num2str(mu),', \sigma=', num2str(sigma)])

% Perform linear regression analysis on growth rate and number of colonies
% versus Nseed actual
%% Plot the expecter variance for each N0 versus the observed variance

% first solve for b+d
for i = 1:length(greensum)
    g = greensum(i).mean_g;
    varg = greensum(i).var_g;
    N0 = greensum(i).Nseed;
    t = green(120).time(end);
    bplusd(i) = ((g*varg)./N0)*(1./(exp(2*g*t)-exp(g*t)));
    d(i) = 0.5*(bplusd(i)-g);
    b(i) = bplusd(i)-d(i);
end








