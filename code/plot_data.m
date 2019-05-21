% process uses mixtures to present and analyze data
close all; clear all; clc
S = load('../out/puref.mat');
pure= S.pure;

%% Further Analysis

% create new structure with g distribution, median, mean, variance, max,
% and percent take-off

% Set a threshold (somewhat arbitrary?) for take-ff
thres = 1000;
for j = 1:length(pure)
    
    if max(pure(j).post_sg) > thres 
        pure(j).takeoff = 1;
    end
    if max(pure(j).post_sg) <thres 
        pure(j).takeoff = 0;
    end
end
% Calculate percent take off for each N plate
ct_takeoff = zeros([1 5]); % 1, 2, 4, 8, 64
ct = zeros(1,5);
for j = 1:length(pure)
    if pure(j).Nplate == 1
        ct(1) = ct(1)+1;
        if pure(j).takeoff == 1
           ct_takeoff(1) = ct_takeoff(1)+1;
        end
    end
    if pure(j).Nplate == 2
        ct(2) = ct(2)+1;
        if pure(j).takeoff == 1
           ct_takeoff(2) = ct_takeoff(2)+1;
        end
    end
    if pure(j).Nplate == 4
        ct(3) = ct(3)+1;
        if pure(j).takeoff == 1
           ct_takeoff(3) = ct_takeoff(3)+1;
        end
    end
    if pure(j).Nplate == 8
        ct(4) = ct(4)+1;
        if pure(j).takeoff == 1
           ct_takeoff(4) = ct_takeoff(4)+1;
        end
    end
    if pure(j).Nplate == 64
        ct(5) = ct(5)+1;
        if pure(j).takeoff == 1
           ct_takeoff(5) = ct_takeoff(5)+1;
        end
    end
end

pct_takeoff = 100*(ct_takeoff./ct); 

% Need to find simple way to come up with summary statistics for g
% distributions in the pure structure
% Convert pure t a table with just N plate and g, then rearrange
% accordingly? Use summarize and then reinput into new structure with g
% distribution
for i = 1:5
    puresum(i).gdistrib = [];
    puresum(i).tkfdistrib = [];
end
for j = 1:length(pure)
    if pure(j).Nplate ==1
        puresum(1).Nplate = 1;
        puresum(1).time = pure(j).time;
        puresum(1).cc = pure(j).cc;
        puresum(1).gdistrib = vertcat(puresum(1).gdistrib, pure(j).lsq_g);
        puresum(1).tkfdistrib = vertcat(puresum(1).tkfdistrib, pure(j).takeoff);
    end
    if pure(j).Nplate ==2
        puresum(2). Nplate = 2;
        puresum(2).time = pure(j).time;
        puresum(2).cc = pure(j).cc;
        puresum(2).gdistrib = vertcat(puresum(2).gdistrib, pure(j).lsq_g);
        puresum(2).tkfdistrib = vertcat(puresum(2).tkfdistrib, pure(j).takeoff);
    end
    if pure(j).Nplate ==4
        puresum(3).Nplate = 4;
        puresum(3).time = pure(j).time;
        puresum(3).cc = pure(j).cc;
        puresum(3).gdistrib = vertcat(puresum(3).gdistrib, pure(j).lsq_g);
        puresum(3).tkfdistrib = vertcat(puresum(3).tkfdistrib, pure(j).takeoff);
    end
    if pure(j).Nplate ==8
        puresum(4). Nplate = 8;
        puresum(4).time = pure(j).time;
        puresum(4).cc = pure(j).cc;
        puresum(4).gdistrib = vertcat(puresum(4).gdistrib, pure(j).lsq_g);
        puresum(4).tkfdistrib = vertcat(puresum(4).tkfdistrib, pure(j).takeoff);
    end
    if pure(j).Nplate ==64
        puresum(5).Nplate = 64;
        puresum(5).time = pure(j).time;
        puresum(5).cc = pure(j).cc;
        puresum(5).gdistrib = vertcat(puresum(5).gdistrib, pure(j).lsq_g);
        puresum(5).tkfdistrib = vertcat(puresum(5).tkfdistrib, pure(j).takeoff);
    end
end

for i = 1:length(puresum)
    puresum(i).pct_takeoff = pct_takeoff(i);
    puresum(i).var_g = var(puresum(i).gdistrib);
    puresum(i).std_g = std(puresum(i).gdistrib);
    puresum(i).mean_g = mean(puresum(i).gdistrib);
   ind = puresum(i).tkfdistrib ==1;
    puresum(i).mean_gtkoff = mean(puresum(i).gdistrib(ind));
    puresum(i).std_gtkoff = std(puresum(i).gdistrib(ind));
    puresum(i).CI_g = [puresum(i).mean_g-2*puresum(i).std_g puresum(i).mean_g+2*puresum(i).std_g];
    puresum(i).model_gavg = singleexpmodel(puresum(i).mean_g,puresum(i).Nplate, puresum(i).time);
    puresum(i).model_gtkoff = singleexpmodel(puresum(i).mean_gtkoff,puresum(i).Nplate, puresum(i).time); 
end
%% Log Nt/No vs time for all cultures and takeoffs
figure;
hold off
for i = 1:length(puresum)
    hold on
    plot(puresum(i).time, log(puresum(i).model_gavg./(puresum(i).Nplate)), 'color',puresum(i).cc, 'LineWidth', 3)
    text(puresum(i).time(end-20), log(puresum(i).model_gavg(end)./(puresum(i).Nplate)),['N_{0}= ', num2str(puresum(i).Nplate), ', g_{mean} = ', num2str(round(puresum(i).mean_g, 4)), ' +/- ', num2str(round(2*puresum(i).std_g, 4))], 'HorizontalAlignment','left','VerticalAlignment','top','color','k')
    xlabel('Time (hrs)')
    ylabel('log(Nt/N_{0})')
    xlim([0 620])
    title ('log(Nt/N_{0}) for all cultures')
end

figure;
hold off
for i = 1:length(puresum)
    hold on
    plot(puresum(i).time, log(puresum(i).model_gtkoff./(puresum(i).Nplate)), 'color',puresum(i).cc, 'LineWidth', 3)
    text(puresum(i).time(end-20), log(puresum(i).model_gtkoff(end)./(puresum(i).Nplate)),['N_{0}= ', num2str(puresum(i).Nplate), ', g_{takeoff} = ', num2str(round(puresum(i).mean_gtkoff, 4)), ' +/- ', num2str(round(2*puresum(i).std_gtkoff, 4))], 'HorizontalAlignment','left','VerticalAlignment','top','color','k')
    xlabel('Time (hrs)')
    ylabel('log(Nt/N_{0})')
    xlim([0 620])
    title ('log(Nt/N_{0}) for takeoffs')
end

%% Per capita growth rate per cell number N colored by initial plating
figure;
hold off
for i = 1:length(pure)
%     if pure(i).Nplate == 4 || pure(i).Nplate == 64
    time = pure(i).time;
    t = round(time,0);
    days = 0:48:t(end);% length of time for calculating per capita g
    ind = find(ismember(t, days));
    hold on
    plot(pure(i).post_sg(ind(1:end-1)), pure(i).percapitag, '.', 'LineWidth', 0.3, 'color', pure(i).cc)
%     end
    xlabel('Cell Number N')
    ylabel('Per Capita Growth Rate')
    ylim([-2, 4])
    xlim([0, 2e2])
end

%% Variance in growth rate normalized

x = [ 1 2 4 8 64];
for i = 1:length(puresum)-1
    y(i) = puresum(i).var_g;
    normy(i) = (x(i).*puresum(i).var_g)/x(5);
end
c = categorical({'N_{0}=1','N_{0}=2','N_{0}=4','N_{0}=8','N_{0}=64'});
figure;
hold off
bar(y)
title ('Variance in growth rate')

hold on
figure;
bar(normy)
title('Normalized variance in growth rate')




%%
figure;
hold off;
for i = 1:length(puresum)
    hold on
    plot(puresum(i).Nplate, puresum(i).pct_takeoff, 'o', 'LineWidth', 5, 'color', puresum(i).cc);
    text(puresum(i).Nplate, puresum(i).pct_takeoff, ['N_{0}= ', num2str(puresum(i).Nplate), ', take-off_{obs} = ', num2str(round(puresum(i).pct_takeoff, 0)), '%'], 'HorizontalAlignment','left','VerticalAlignment','top','color','k')
    xlabel ('N0 plated')
    ylabel('Percent Takeoff')
    title('Percent Takeoff vs. N0')
    xlim([0 80])
    ylim ([0 100])
end
hold on
for i = 1:length(puresum)
    expected(i,:)=[puresum(i).Nplate puresum(1).pct_takeoff*puresum(i).Nplate];
    plot(puresum(i).Nplate, puresum(1).pct_takeoff*puresum(i).Nplate, '*', 'color', puresum(i).cc)
    text(puresum(i).Nplate, puresum(1).pct_takeoff*puresum(i).Nplate, ['N_{0}= ', num2str(puresum(i).Nplate),' take-off_{exp} = ', num2str(round(puresum(1).pct_takeoff*puresum(i).Nplate, 0)), '%'], 'HorizontalAlignment','left','VerticalAlignment','top','color',puresum(i).cc)
end 
plot(expected(:,1), expected(:,2), '-', 'LineWidth', .2, 'color', 'k')

text(puresum(i).Nplate, puresum(1).pct_takeoff*puresum(i).Nplate, ['take-off_{exp} = ', num2str(round(puresum(1).pct_takeoff*puresum(i).Nplate, 0)), '%'], 'HorizontalAlignment','left','VerticalAlignment','top','color','k')

save('../out/purefp.mat', 'pure')
save('../out/puresum.mat', 'puresum')

%% Plot results of various smoothing algorithms,


figure;
hold off
start = 1;
last = length(pure);
subplot(1,4,1)
for i = start:last%:length(pure)
    cc = pure(i).cc;
    semilogy(pure(i).time, pure(i).cellnum, 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Spaghetti Plot of Raw Data')
    xlim([0 504])
    ylim([1 10e5])

end

subplot(1,4,2)
for i = start:last%1:length(pure)
    cc = pure(i).cc;
    semilogy(pure(i).time, pure(i).post_sg, 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title(' Savitsky-Golay Filter, span =30 ')
    xlim([0 504])
    ylim([1 10e5])
   
end

subplot(1,4,3)
for i = start:last%1:length(pure)
    cc = pure(i).cc;
    semilogy(pure(i).time, (pure(i).post_rl), 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title(' Local Regression')
    xlim([0 504])
    ylim([1 10e5])
   
end
subplot(1,4,4)
for i = start:last%1:length(pure)
    cc = pure(i).cc;
    semilogy(pure(i).time, (pure(i).post_rl10), 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Local Regression , span =10 ')
    xlim([0 504])
    ylim([1 10e5])
   
end
%% Just look at Smith- Golay and 8 vs 1
start = 1
last = length(pure)
figure; hold off
for i = start:last%1:length(pure)
    if pure(i).Nplate ==1
        cc = pure(i).cc;
        semilogy(pure(i).time, pure(i).post_sg, 'color', cc)
    end
    if pure(i).Nplate == 4
        cc = pure(i).cc;
        semilogy(pure(i).time, pure(i).post_sg, 'color', cc)
    end
%     if pure(i).Nplate ==1
%         cc = pure(i).cc;
%         plot(pure(i).time, pure(i).post_sg, 'color', 'r')
%     end
%     if pure(i).Nplate == 64
%         cc = pure(i).cc;
%         plot(pure(i).time, pure(i).post_sg, 'color', 'b')
%     end
        
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Smoothed Data: Savitsky-Golay Filter, span =30 ')
    xlim([0 504])
     ylim([1 10e4])
   
end
%% Plot single exponential model (Bayesian and LSQnonlin) vs. smoothed data
figure;
hold off
for j = 200:210%:length(pure)
        pause
        plot(pure(j).time, pure(j).post_sg, 'k*')
        hold on
        plot(pure(j).time, pure(j).Nmodel, 'g-', 'LineWidth', 2)
        plot(pure(j).time, pure(j).Nmodellsq, 'b-', 'LineWidth', 2)
        hold off
 
end
%% Plot all on one graph
% first do lsqnonlin fits

figure; hold off
start = 1;
last = length(pure);
for i = start:last%1:length(pure)
    if pure(i).Nplate ==1
        cc = pure(i).cc;
        plot(pure(i).time, log(pure(i).Nmodellsq./pure(i).Nplate), 'color', 'r')
    end
    if pure(i).Nplate == 2
        cc = pure(i).cc;
%         plot(pure(i).time, log(pure(i).Nmodellsq./pure(i).Nplate), 'color', 'b')
    end
    if pure(i).Nplate ==4
        cc = pure(i).cc;
%         plot(pure(i).time, log(pure(i).Nmodellsq./pure(i).Nplate), 'color', 'g')
    end
    if pure(i).Nplate == 8
        cc = pure(i).cc;
        plot(pure(i).time, log(pure(i).Nmodellsq./pure(i).Nplate), 'color', 'm')
    end
    
     if pure(i).Nplate == 64
        cc = pure(i).cc;
        plot(pure(i).time, log(pure(i).Nmodellsq./pure(i).Nplate), 'color', 'k')
    end
        
    hold on
    xlabel ('time (hours)')
    ylabel('Log(Nt/No)')
    title('Model fits using lsqnonlin')
    xlim([0 504])
%     ylim([0 3000])
   
end
%% Using Bayesian
figure; hold off
start = 1;
last = length(pure);
for i = start:last%1:length(pure)
    hold on
    if pure(i).Nplate ==1
        cc = pure(i).cc;
        plot(pure(i).time, (pure(i).Nmodel), 'color', 'r')
    end
    if pure(i).Nplate == 2
        cc = pure(i).cc;
        plot(pure(i).time, (pure(i).Nmodel), 'color', 'b')
    end
    if pure(i).Nplate ==4
        cc = pure(i).cc;
        plot(pure(i).time, (pure(i).Nmodel), 'color', 'g')
    end
    if pure(i).Nplate == 8
        cc = pure(i).cc;
        plot(pure(i).time, (pure(i).Nmodel), 'color', 'm')
    end
    
     if pure(i).Nplate == 64
        cc = pure(i).cc;
        plot(pure(i).time, log(pure(i).Nmodel), 'color', 'k')
    end
        
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Model fits using Bayesian')
    xlim([0 504])
%     ylim([0 3000])
   
end

%% Growth rate as a function of number plates
figure;
hold off;
for i= 1:length(pure)
    hold on
    if pure(i).Nplate ==1
        cc = pure(i).cc;
        plot(pure(i).Nplate, pure(i).lsq_g,  'r*')
    end
    if pure(i).Nplate == 2
        cc = pure(i).cc;
        plot(pure(i). Nplate, pure(i).lsq_g, 'b*')
    end
    if pure(i).Nplate ==4
        cc = pure(i).cc;
        plot(pure(i).Nplate, pure(i).lsq_g, 'g*')
    end
    if pure(i).Nplate == 8
        cc = pure(i).cc;
        plot(pure(i).Nplate, pure(i).lsq_g, 'm*')
    end
     if pure(i).Nplate == 64
        cc = pure(i).cc;
        plot(pure(i).Nplate, pure(i).lsq_g, 'k*')
     end
     hold on
     
     xlabel ('N0')
     ylabel ('growth rate g')
     title ('Growth rate as a function of N0')
end
