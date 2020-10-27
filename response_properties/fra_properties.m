function [tc, fraproperties, psth, peakdelay, rpw, n0] = fra_properties(spk, raster, trigger, params)
% Inputs:
%   spk: units from spike sorting with spiketimes in the structure
%   raster: raster of spikes in response to pure tones
%
% Outputs:
%   tc: tuning curve, calculated with significant timebins (> rest + 3std)
%   fraproperties: boundry of FRA, threshold, Q10 etc.
%       boundry of FRA is generated from FRA smoothed by 3x3 Gaussian
%       filter, baseline as last 100ms mean + 2std (without smoothing)
%   psth: origional psth from rastr without smooth
%   delay: first timebin with significant increase, so does the following timebin
%       from smoothed psth
%   rpw: response window, time interval with significant increase in psth
%       from smoothed psth
%   n0: number of total spikes (first to last trigger + 300ms)
%   
%   first, the significance of psth is determined by a moving 5 ms wondow
%   fr5ms > mean(psth(201:300)) + 3 * std(psth(201:300))
%   for significant psth, psth is smoothed and continuous significant time
%   bins around the peak is used as response window (rpw)
% 
% First, the significance of PSTH is determined by signed rank test of the
% number of spikes during sound and before sound
% 
% adapted from FRA_properties_NH, 12/17/2019
% Congcong, 12/20/2019
% set variables
levels = unique(raster.atten); % row1 = 70 dB SPL
freqs = unique(raster.freq);
nlevels = length(levels);
nfreqs = length(freqs);
trange = [0 300];
binsize = 1;
allTimes = [];
edges = trange(1):binsize:trange(2);
lim = [0 50];

% ----------------to determine the significance o psth (signed-rank test, p=0.01)--------------------
spkcount = zeros(size(raster.rastermat,1),2);
for ii = 1:size(raster.rastermat,1)
    eachTimes = raster.rastermat{ii};
    eachTimes = eachTimes(:);%spike time to each tone (stimulus presentation + silence)
    allTimes = [allTimes; eachTimes];
    
    f =find(eachTimes>= lim(1) & eachTimes < lim(2));
    spkcount(ii,1) = length(f); % number of rasterspikes during stimulus presentation
    f =find(eachTimes>= edges(end)-lim(2) & eachTimes < edges(end) - lim(1));
    spkcount(ii,2) = length(f); % number of spikes during silence
end

[~,significance] = signrank(spkcount(:,1),spkcount(:,2),'alpha', 0.01); %,'tail','right');

[psth,~] = histcounts(allTimes,edges);
n0 = length(allTimes); %number of spikes

% --------------------get significant window-------------------------------
% if psth is significant, get time window with significant increase in firing rate
% smooth setting: 5 bins average
if significance
    psth_smooth = smooth(psth, 5);
    rest = mean(psth_smooth(201:300));
    sigma = std(psth_smooth(201:300));
    sigbins = find(psth_smooth > rest + 3 * sigma);
    if ~isempty(sigbins)
        [~, peak_idx] = max(psth_smooth); %index of peak in psth
        peakdelay = peak_idx;
        peak_idx = find(sigbins == peak_idx); %index of peak in sigbins
        ii = peak_idx;
        while ii > 1
            if sigbins(ii ) - sigbins(ii - 1) == 1
                ii = ii - 1;
            else
                break
            end
        end
        rpw(1) = sigbins(ii)-1; %start of significant increase in firing rate in two continuous time bins
        ii = peak_idx;
        while ii < length(sigbins)
            if sigbins(ii + 1) - sigbins(ii) == 1
                ii = ii+1;
            else
                break
            end
        end
        rpw(2) = sigbins(ii); %end of significant increase in firing rate in two continuous time bins
    else
        peakdelay = NaN;
        rpw = [0 50];
    end
else
    peakdelay = NaN;
    rpw = [0 50];
end

%calculate tuning curve and decide tc boundary
tc = calculate_tuning_curve(spk,trigger,params,rpw);
ntc = calculate_tuning_curve(spk,trigger,params,[300-diff(rpw), 300]); %get null fra with spontaneous activity
ntc = ntc.tcmat(:);
base = mean(ntc) + 3*std(ntc);

spikes = tc.tcmat;
spikes(spikes < base) = 0;
spikes = smoothmat(spikes, 2);

bounds=(nlevels+1)*ones(nfreqs, 1);
%calculate a threshold based on null tc calculated with spontaneous activity.
for ii=1:nfreqs
    for jj=1:nlevels
        if spikes(jj, ii)
            if (jj<nlevels-3 && spikes(jj+1, ii) && spikes(jj+2, ii) && spikes(jj+3, ii) && spikes(jj+4, ii))||...
               (jj == nlevels-3 && spikes(jj+1, ii) && spikes(jj+2, ii) && spikes(jj+3, ii)) ||...
               (jj == nlevels-2 && spikes(jj+1, ii) && spikes(jj+2, ii) )||...
               (jj == nlevels-1 && spikes(jj+1, ii))
                if ii == 1 || spikes(jj, ii-1)>0
                    bounds(ii) = jj + 1;
                    break
                end
                 
            end
        end
        
    end
end
for ii=2:length(bounds)-1 % remove odd peaks
    if bounds(ii-1)==nlevels+1 && bounds(ii+1)==nlevels+1 && bounds(ii)<nlevels-5
        bounds(ii)=nlevels+1;
    end
end

% Call FRApropertiesCx01
% extract fraproperties (original script is FRApropertiesTH05 from DPhil project)
%spikes = smoothmat(tc.tcmat, 2);
spikes = tc.tcmat;
spikes(spikes <= base) = 0;

[frapropertiestmp]=FRApropertiesCx01(bounds,freqs,spikes,levels); %,fileName,channel,idx_window)

% add outputs to fraproperties
fraproperties.threshold = frapropertiestmp.threshold;
fraproperties.BF = frapropertiestmp.BF;
fraproperties.CF = frapropertiestmp.CF;
fraproperties.BW10 = frapropertiestmp.BW10;
fraproperties.Q10 = frapropertiestmp.Q10;
fraproperties.BW30 = frapropertiestmp.BW30;
fraproperties.Q30 = frapropertiestmp.Q30;
fraproperties.BW40 = frapropertiestmp.BW40;
fraproperties.Q40 = frapropertiestmp.Q40;
fraproperties.boundry = bounds; %single peak boundry containing bf
fraproperties.sboundry = frapropertiestmp.newbounds; %single peak boundry containing bf
fraproperties.base = base;
end







