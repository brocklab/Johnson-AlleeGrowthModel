% process uses mixtures to present and analyze data
close all; clear all; clc
S = load('../out/pure.mat');
pure= S.pure;


%% Plot results of various smoothing algorithms,


figure;
hold off
start = 1;
last = length(pure);
subplot(1,4,1)
for i = start:last%:length(pure)
    cc = pure(i).cc;
    plot(pure(i).time, pure(i).cellnum, 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Spaghetti Plot of Raw Data')
    xlim([0 550])
    ylim([0 3000])

end

subplot(1,4,2)
for i = start:last%1:length(pure)
    cc = pure(i).cc;
    plot(pure(i).time, pure(i).post_sg, 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Smoothed Data: Savitsky-Golay Filter, span =30 ')
    xlim([0 550])
    ylim([0 3000])
   
end

subplot(1,4,3)
for i = start:last%1:length(pure)
    cc = pure(i).cc;
    plot(pure(i).time, pure(i).post_rl, 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Smoothed Data: Local Regression')
    xlim([0 550])
    ylim([0 3000])
   
end
subplot(1,4,4)
for i = start:last%1:length(pure)
    cc = pure(i).cc;
    plot(pure(i).time, pure(i).post_rl10, 'color', cc)
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Smoothed Data: Local Regression , span =10 ')
    xlim([0 550])
    ylim([0 3000])
   
end
%% Just look at Smith- Golay and 8 vs 1
start = 1
last = length(pure)
figure; hold off
for i = start:last%1:length(pure)
    if pure(i).Nplate ==1
        cc = pure(i).cc;
        plot(pure(i).time, pure(i).post_sg, 'color', cc)
    end
    if pure(i).Nplate == 8
        cc = pure(i).cc;
        plot(pure(i).time, pure(i).post_sg, 'color', 'g')
    end
    if pure(i).Nplate ==1
        cc = pure(i).cc;
        plot(pure(i).time, pure(i).post_sg, 'color', 'r')
    end
    if pure(i).Nplate == 64
        cc = pure(i).cc;
        plot(pure(i).time, pure(i).post_sg, 'color', 'b')
    end
        
    hold on
    xlabel ('time (hours)')
    ylabel('Cell Number')
    title('Smoothed Data: Savitsky-Golay Filter, span =30 ')
    xlim([0 504])
%     ylim([0 3000])
   
end


