%PLOT_PSTH plot ratser data as psth
% PLOT_PSTH(RASTER) draw raster plot stored in the raster struct (see calculate_raster for fields).
% Natsumi 240217

function [significance,p] = plot_psth(raster,lim)

if nargin < 2
    lim = [0 50];
end

%Checking input data
if(length(raster) > 1)
    error('Can only plot one tuning curve at a time');
end

binsize = 5; % 2 ms

%Plotting psth
hold on;
allTimes = [];
edges = 0:binsize:raster.Dtrigger+1;

for i = 1:size(raster.rastermat,1)
    eachTimes = raster.rastermat{i};
    allTimes = [allTimes; eachTimes];
    
    f =find(eachTimes>= lim(1) & eachTimes < lim(2));
    spkcount(i,1) = length(f); % number of spikes during stimulus presentation
    f =find(eachTimes>= edges(end)-lim(2) & eachTimes < edges(end) - lim(1));
    spkcount(i,2) = length(f); % number of spikes during silence
    
end
[p,significance] = signrank(spkcount(:,1),spkcount(:,2),'alpha', 0.001); %,'tail','right');
[N,edges] = histcounts(allTimes,edges);
h_hist = bar(edges(1:end-1)+edges(2)/2, N, 'k');
y = get(gca, 'ylim'); y = y(2);
xlim([0 edges(end)])
set(h_hist, 'barwidth', 1)
% ylabel('Spike count');
% xlabel('Time (ms)');
end