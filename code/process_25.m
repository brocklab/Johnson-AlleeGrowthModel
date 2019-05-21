% Process individual well data and wells grouped by Nseed data

% load fitted and processed data
close all; clear all; clc
S = load('../out/greencf.mat');
greenc= S.greenc;
S = load('../out/greencsum.mat');
greencsum= S.greencsum;
%%
% plot raw data by color
figure (1)
subplot(2,2,1:2)
for j = 1:length(greenc)
    if greenc(j).Nseed == 2
    plot(greenc(j).time, greenc(j).cellnum, 'color', 'r', 'LineWidth', 1)
    end
    if greenc(j).Nseed == 5
    plot(greenc(j).time, greenc(j).cellnum, 'color', 'g', 'LineWidth', 1)
    end
    hold on
end
xlim([0, greenc(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('All wells')
subplot(2,2,3)
for j = 1:length(greenc)
    if greenc(j).Nseed ==2
    plot(greenc(j).time, greenc(j).cellnum, 'color', 'r', 'LineWidth', 1)
    hold on
    end
end
xlim([0, greenc(1).time(end)])
ylim([0, 1000])
xlabel('Time (hours)')
ylabel('N')
title(['N_{seed} = 2, Pct_{take-off}= ', num2str(round(greencsum(1).pct_takeoff,0)), '%'])

subplot(2,2,4)
for j = 1:length(greenc)
    if greenc(j).Nseed ==5
    plot(greenc(j).time, greenc(j).cellnum, 'color', 'g', 'LineWidth', 1)
    hold on
    end
end
xlim([0, greenc(1).time(end)])
xlabel('Time (hours)')
ylabel('N')
title(['N_{seed} = 5, Pct_{take-off}= ', num2str(round(greencsum(2).pct_takeoff,0)), '%'])
%%

figure(2)
hold off
subplot(1,3,1)
for j = 1:length(greenc)
    time = greenc(j).time;
    if greenc(j).takeoff == 1
        if greenc(j).Nseed == 2
        plot(greenc(j).time, log(greenc(j).cellnum/greenc(j).Nseed), 'color', 'r', 'LineWidth', 1)
        end
        if greenc(j).Nseed == 5
        plot(greenc(j).time, log(greenc(j).cellnum/greenc(j).Nseed), 'color', 'g', 'LineWidth', 1)
        end
    end
    hold on
end
xlim([0, greenc(1).time(end)])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Raw data takeoffs only')


subplot(1,3,2)
for j = 1:length(greenc)
    time = greenc(j).time;
    if greenc(j).takeoff == 1
        if greenc(j).Nseed == 2
        plot(greenc(j).time, log(greenc(j).Nmodellsq/greenc(j).Nseed), 'color', 'r', 'LineWidth', 1)
        end
        if greenc(j).Nseed == 5
        plot(greenc(j).time, log(greenc(j).Nmodellsq/greenc(j).Nseed), 'color', 'g', 'LineWidth', 1)
        end
    end
    hold on
end
xlim([0, greenc(1).time(end)])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Per Capita Growth takeoffs')



subplot(1,3,2)
hold off
for j = 1:length(greencsum)
    model_low =  singleexpmodel(greencsum(j).mean_gtkoff - 1.96*sqrt(greencsum(j).var_gtkoff),greencsum(j).Nseed, time);
    model_high =  singleexpmodel(greencsum(j).mean_gtkoff + 1.96*sqrt(greencsum(j).var_gtkoff),greencsum(j).Nseed, time);
    if greencsum(j).Nseed == 2
        plot(time, greencsum(j).model_gtkoff, 'color', 'r', 'LineWidth', 3)
    end
    if greencsum(j).Nseed == 5
        plot(time, greencsum(j).model_gtkoff, 'color', 'g', 'LineWidth', 3)
    end
    hold on
    %plot(time, model_low, 'color', greensum(j).cc, 'LineWidth', 1)
    %plot(time, model_high, 'color', greensum(j).cc, 'LineWidth', 1)
   % text(time(end), log(greensum(j).model_gtkoff(end)./greensum(j).Nseed), ['N_{0}= ', num2str(greensum(j).Nseed), ', g_{mean}= ', num2str(round(greensum(j).mean_gtkoff,4))], 'HorizontalAlignment', 'Right')
    hold on
end

xlim([0, greencsum(j).time(end)])
xlabel('Time (hours)')
ylabel('N')
title('Mean population growth dynamics')
tint = [99 198 300];
for i = 1:3
ind(i) = find(ismember( round(time,0), tint(i)));
end

subplot(1,3,3)
for j = 1:length(greencsum)
    time = greencsum(j).time;
    model_low =  singleexpmodel(greencsum(j).mean_gtkoff - sqrt(greencsum(j).var_gtkoff),greencsum(j).Nseed, time);
    model_high =  singleexpmodel(greencsum(j).mean_gtkoff + sqrt(greencsum(j).var_gtkoff),greencsum(j).Nseed, time);
    
    if greencsum(j).Nseed ==2
    plot(time, log(greencsum(j).model_gtkoff./greencsum(j).Nseed), 'color', 'r', 'LineWidth', 3)
    hold on
    errorbar([100;200;300], log(greencsum(j).model_gtkoff(ind)./greencsum(j).Nseed), 0.5.*(log(model_high(ind)./greencsum(j).Nseed)-log(model_low(ind)./greencsum(j).Nseed)), 'color', 'r', 'LineWidth', 1.3)
    end
    
    if greencsum(j).Nseed ==5 
    plot(time, log(greencsum(j).model_gtkoff./greencsum(j).Nseed), 'color', 'g', 'LineWidth', 3)
    hold on
    errorbar([100;200;300], log(greencsum(j).model_gtkoff(ind)./greencsum(j).Nseed), 0.5.*(log(model_high(ind)./greencsum(j).Nseed)-log(model_low(ind)./greencsum(j).Nseed)), 'color', 'g', 'LineWidth', 1.3)
    end
    %plot(time, log(model_low/greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 1)
    %plot(time, log(model_high/greensum(j).Nseed), 'color', greensum(j).cc, 'LineWidth', 1)
    %text(time(end), log(greencsum(j).model_gtkoff(end)./greencsum(j).Nseed), ['N_{0}= ', num2str(greencsum(j).Nseed), ', g_{mean}= ', num2str(round(greencsum(j).mean_gtkoff,4)),' +/- ', num2str(round(sqrt(greencsum(j).var_gtkoff),3)),' cells/hr'], 'HorizontalAlignment', 'Right')
    hold on
end

xlim([0, greenc(1).time(end)])
xlabel('Time (hours)')
ylabel('log(N/N_{0})')
title('Average Per Capita Growth by N_{seed}')
%%
for j = 1:length(greensum)
    norm_var(j) = greensum(j).norm_var
end
%%
figure;
bar(norm_var)
ylabel('Normalized Variance')
%%
x1 = ones(length(greencsum(1).gtkoff));
x2 = 2*ones(length(greencsum(2).gtkoff));
g1 = greencsum(1).gtkoff;
g2 = greencsum(2).gtkoff;

gr1 = greencsum(1).gtot;
gr2 = greencsum(2).gtot;

figure;
subplot(1,2,1)
plot(x1, g1, 'ro')
hold on
plot(x2, g2, 'go')
xlim([ 0 3])
title('Take-off growth rates for 2 & 5 cells')

xr1 = ones(length(greencsum(1).gtot));
xr2 = 2*ones(length(greencsum(2).gtot));
subplot(1,2,2)
plot(xr1, gr1, 'ro')
hold on
plot(xr2, gr2, 'go')
xlim([ 0 3])
title('All growth rates for 2&5 cells')
%%
figure;




    


