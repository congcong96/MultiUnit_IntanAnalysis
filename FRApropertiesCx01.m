function [fraproperties]=FRApropertiesCx01(bounds,freqs,spikes,levels)
%% FRAproperties05
% this script is for calculating threshold, BF, CF, Q10, Q40, ratio of tuning curve (Tratio) and peaks of tuning curve (NumPeak)
% woking with following the script drawFRAbounds03
% original script is FRAproperties01
% ratio of tuning curve (Tratio) - see Morel et al., 1987 Tratio NumPeak

% Oct 2013 Natsumi

% changed FRAproperties05 for AnalyzeRecordingDataTH_01.m
% Oct 2014 Natsumi

% changed FRApropertiesTH01 some lines were taken from original
% FRAproperties (e.g. bw bandwidth bandwidth40)
% this was checked in 
% August 2015 Natsumi

% changed FRApropertiesTH03 the same as FRApropertiesTH04 from FRApropertiesTH02.
% get newbounds
% March 2017 Natsumi

% updated FRApropertiesTH05 for intan system
% Natsumi 23April2018

%% extract BF (best frequency) - a frequency with the maximum response
% nfreqs=size(unique(freqs),1);
% freqspace=linspace(min(log(freqs)),max(log(freqs)),nfreqs*2);
% maxspikes=max(max(spikes));
% [f1 f2]=find(spikes==maxspikes);
% BFvals=mean(freqspace(f2*2));
% BF=exp(BFvals)/1000;

%% extract BF2 (another dbest frequency) - a frequency with the maximum response across all levels
nfreqs=size(unique(freqs),1);
freqspace=linspace(min(log(freqs)),max(log(freqs)),nfreqs*2);
acrosslevels=sum(spikes,1);
[~, ff]=max(acrosslevels);
BFvals2=mean(freqspace(ff*2));
BF2=exp(BFvals2)/1000;

%%  if there are frequency gaps in the lowest intensity, take a part that has BF.
[I2]=find(bounds<(length(unique(levels))+1)); % find the bound in FRA
checkI2 = diff(I2);
fI2 = find(checkI2>1);
if ~isempty(fI2)
    newbounds = ones(size(bounds))*(length(unique(levels))+1);
    for i = length(fI2)+1:-1:1
        if i == 1
            startidx(i) = I2(1);
            endidx(i) = I2(fI2(i));
        elseif i == length(fI2)+1
            startidx(i) = I2(fI2(i-1)+1);
            endidx(i) =  I2(end);
        else
            startidx(i) = I2(fI2(i-1)+1);
            endidx(i) =  I2(fI2(i));
        end
    end
    
    
    for i = length(startidx)
        if startidx(i) <= mean(ff)  && endidx(i) >= mean(ff)
            newbounds(startidx(i):endidx(i)) = bounds(startidx(i):endidx(i));
        end
    end
    if exist('newbounds','var')
        % calculate area
        for i = length(startidx):-1:1
            areas(i) = sum((length(unique(levels))+1)-bounds(startidx(i):endidx(i)));
        end
        maxareaidx = find(areas == max(areas));
        newbounds(startidx(maxareaidx):endidx(maxareaidx)) = bounds(startidx(maxareaidx):endidx(maxareaidx));
    end
else
    newbounds = bounds;
end

%% Extract threshold as the minimum intensity to evoke a response
if length(unique(newbounds))~=1
    threshold = levels(min(newbounds));
else
    threshold = NaN;
end

%% Extract Q10
% Q10 = CF / the width of the tuning curve 10 dB above threshold 

q=(min(newbounds)+2);
[I]=find(newbounds<=q);
    
index=1;
ind=1;
bw=[];
tuning=0;

if isempty(I) || min(newbounds)==(length(unique(levels))+1)% || isnan(CF) %isempty(I) || min(newbounds)==5 || isnan(CF)
    Q10=NaN;
    bandwidth=NaN;
    Tratio=NaN;
    tuning=NaN;
    highfreq=NaN;
    lowfreq=NaN;
    Tratio=NaN;
    
else
    if length(I)==1
        bw=[I I];
    else
        for ii=1:length(I)
            if ii==length(I)
                if I(end)-1~=I(end-1)
                    bw=[bw;I(end) I(end)];
                elseif I(end)-1==I(end-1)
                    bw=[bw;I(index) I(end)];
                end
                
            elseif I(ii+1)~=I(ii)+1
                bw=[bw; I(index) I(ii)];
                index=ii+1;
                ind=ind+1;
            end
        end
    end
    
    bandwidth=[];
    for ii=1:size(bw,1)
        
        if bw(ii,1)==bw(ii,2)
            if bw(ii,1)==nfreqs
                bandwidth=[bandwidth;exp(freqspace(bw(ii,1)*2))/1000-exp(freqspace(bw(ii,1)*2-2))/1000];
            else
                bandwidth=[bandwidth;exp(freqspace(bw(ii,1+1)*2))/1000-exp(freqspace(bw(ii,1)*2-1))/1000];
            end
        else
            bandwidth=[bandwidth;exp(freqspace(bw(ii,2)*2))/1000-exp(freqspace(bw(ii,1)*2))/1000];
        end
    end

    highfreq=exp(freqspace(bw(size(bw,1),2)*2))/1000;
    lowfreq=exp(freqspace(bw(1,1)*2))/1000;
    Tratio=highfreq/lowfreq;
end

BW10freq(1)=highfreq;
BW10freq(2)=lowfreq;
bandwidth=nansum(bandwidth);
BW10=bandwidth;
%Q10=CF/bandwidth;

%% extract CF (characteristic frequency) - a frequency with the lowest intensity to evoke a response
% if there are frequency gaps in the lowest intensity, take a part that has BF.

[I]=find(newbounds==(min(newbounds))); % find the freqs with the lowest intensity.
%freqspace=linspace(min(log(freqs)),max(log(freqs)),nfreqs*2);
gap=[];
idx_gap=1;

if length(I)==1
    gap=[I I];
else
    for ig=1:length(I)
        if ig==length(I)
            if I(end)-1~=I(end-1)
                gap=[gap;I(end) I(end)];
            elseif I(end)-1==I(end-1)
                gap=[gap;I(idx_gap) I(end)];
            end
            
        elseif I(ig+1)~=I(ig)+1
            gap=[gap; I(idx_gap) I(ig)];
            idx_gap=ig+1;
        end
    end
end

numpeak=size(gap,1);

if isempty(I)
    CF=NaN;
elseif length(I)==length(newbounds) && min(newbounds)==(length(unique(levels))+1) %min(newbounds)==5
    CF=NaN;
elseif min(newbounds)==0
    CF=NaN;
else

if size(gap,1)==1
    CFvals=mean(freqspace((gap(1,1):gap(1,end))*2));
    CF=exp(CFvals)/1000;
    CFvals2plot=mean(gap(1,1):gap(1,2));
    %line([mean(gap(1,1):gap(1,2)), mean(gap(1,1):gap(1,2))],[0 (length(unique(levels))+1)],'Color','w'); %line([mean(ff) mean(ff)],[0 5],'Color','w');
else
    % logarithmically weighted mean
    for ii=1:size(gap,1)
        cfvals(ii)=mean(freqspace((gap(ii,1):gap(ii,end))*2));
        numcfvals(ii)=gap(ii,2)-gap(ii,1)+1;
        numcfvalsxcfvals(ii)=cfvals(ii)*numcfvals(ii);
        % for plot
        Pcfvals(ii)=mean((gap(ii,1):gap(ii,end)));
        Pnumcfvals(ii)=gap(ii,2)-gap(ii,1)+1;
        Pnumcfvalsxcfvals(ii)=Pcfvals(ii)*Pnumcfvals(ii);
    end
    
    CFvals=sum(numcfvalsxcfvals)/sum(numcfvals);
    CF=exp(CFvals)/1000;   
    CFvals2plot=sum(Pnumcfvalsxcfvals)/sum(Pnumcfvals);
    %CF2plot=exp(freqspace(round(CFvals2plot*2)))/1000; Thid does not match to CF because freqspace is sparce vector.
    %line([CFvals2plot CFvals2plot],[0 (length(unique(levels))+1)],'Color','w'); %line([mean(ff) mean(ff)],[0 5],'Color','w');

end
end

if ~exist('CFvals2plot')
    CFvals2plot=NaN;
end

% %evaluate Tratio according to the criteria: Morel et al., 1987
%     if CF >= 4
%         if Tratio >= 4
%             tuning=2;%'Broad';
%         elseif Tratio < 4
%             tuning=1;%'Narrow';
%         elseif isnan(Tratio)
%             tuning=0;%'NaN';
%         else
%             disp('error')
%             tuning=0;
%         end
%     elseif CF < 4
%         if Tratio >= 2
%             tuning=2;%'Broad';
%         elseif Tratio < 2
%             tuning=1;%'Narrow';
%         elseif isnan(Tratio)
%             tuning=0;%'NaN';
%         else
%             disp('error')
%             tuning=0;
%         end
%     elseif isnan(CF)
%         tuning=0;%'NaN';
%     end

    if ~isnan(CF) & ~isnan(BW10)
        Q10=CF/BW10;
    else
        Q10=NaN;
    end
    
%% Extract Q30 (repeate for Q30)
% Q30 = CF / the width of the tuning curve 30dB above threshold
q=(min(newbounds)+6); 
[I]=find(newbounds<=q);

index=1;
ind=1;
bw=[];

if isempty(I) || isnan(CF) || min(newbounds)>=16-6 %min(newbounds)>=3
    Q30=NaN;
    bandwidth30=NaN;
    BW30freq(1)=NaN;
    BW30freq(2)=NaN;
else
    if length(I)==1
        bw=[I I];
    else
        for ii=1:length(I)
            if ii==length(I)
                if I(end)-1~=I(end-1)
                    bw=[bw;I(end) I(end)];
                elseif I(end)-1==I(end-1)
                    bw=[bw;I(index) I(end)];
                end
                
            elseif I(ii+1)~=I(ii)+1
                bw=[bw; I(index) I(ii)];
                index=ii+1;
                ind=ind+1;
            end
        end
    end
    
    bandwidth30=[];
    for ii=1:size(bw,1)
        if bw(ii,1)==bw(ii,2)
            if bw(ii,1)==nfreqs
%                 keyboard;
%                 bandwidth40=[bandwidth40;exp(freqspace(bw(ii,1)*2))/1000-exp(freqspace(bw(ii,1)*2-2))/1000];
%             else
                bandwidth30=[bandwidth30;exp(freqspace(bw(ii,2)*2))/1000-exp(freqspace(bw(ii,1)*2-1))/1000];
            end
        else
            bandwidth30=[bandwidth30;exp(freqspace(bw(ii,2)*2))/1000-exp(freqspace(bw(ii,1)*2))/1000];
        end
    end
    bandwidth30=sum(bandwidth30);
    BW30freq(2)=exp(freqspace(bw(size(bw,1),2)*2))/1000;
    BW30freq(1)=exp(freqspace(bw(1,1)*2))/1000;
end
if isempty(bandwidth30)
    bandwidth30=NaN;
end

BW30=bandwidth30;
Q30=CF/BW30;

%% Extract Q40 (repeate for Q40)
% Q40 = CF / the width of the tuning curve 40dB above threshold
% there was no significant differences in Q40 between A1 and VAF. try Q40
% as Polly et al. 2007 used Q42

q=(min(newbounds)+8); 
[I]=find(newbounds<=q);

index=1;
ind=1;
bw=[];

if isempty(I) || isnan(CF) || min(newbounds)>=16-8 %min(newbounds)>=3
    Q40=NaN;
    bandwidth40=NaN;
    BW40freq(1)=NaN;
    BW40freq(2)=NaN;
else
    if length(I)==1
        bw=[I I];
    else
        for ii=1:length(I)
            if ii==length(I)
                if I(end)-1~=I(end-1)
                    bw=[bw;I(end) I(end)];
                elseif I(end)-1==I(end-1)
                    bw=[bw;I(index) I(end)];
                end
                
            elseif I(ii+1)~=I(ii)+1
                bw=[bw; I(index) I(ii)];
                index=ii+1;
                ind=ind+1;
            end
        end
    end
    
    bandwidth40=[];
    for ii=1:size(bw,1)
        if bw(ii,1)==bw(ii,2)
            if bw(ii,1)==nfreqs
%                 keyboard;
%                 bandwidth40=[bandwidth40;exp(freqspace(bw(ii,1)*2))/1000-exp(freqspace(bw(ii,1)*2-2))/1000];
%             else
                bandwidth40=[bandwidth40;exp(freqspace(bw(ii,2)*2))/1000-exp(freqspace(bw(ii,1)*2-1))/1000];
            end
        else
            bandwidth40=[bandwidth40;exp(freqspace(bw(ii,2)*2))/1000-exp(freqspace(bw(ii,1)*2))/1000];
        end
    end
    bandwidth40=sum(bandwidth40);
    BW40freq(2)=exp(freqspace(bw(size(bw,1),2)*2))/1000;
    BW40freq(1)=exp(freqspace(bw(1,1)*2))/1000;
end

if isempty(bandwidth40)
    bandwidth40=NaN;
end

BW40=bandwidth40;
Q40=CF/BW40;

%% organize extracted data
fraproperties.threshold = threshold;
fraproperties.BF = BF2;
fraproperties.CF = CF;
fraproperties.BW10 = BW10freq;
fraproperties.Q10 = Q10;
fraproperties.BW30 = BW30freq;
fraproperties.Q30 = Q30;
fraproperties.BW40 = BW40freq;
fraproperties.Q40 = Q40;
%fraproperties.Tratio = Tratio;
%fraproperties.numpeak = numpeak;
%fraproperties.tuning = tuning;

% for plotting
% fraproperties.CFvals2plot = CFvals2plot;
% fraproperties.BW102plot = BW10;
% fraproperties.BW402plot = BW40;
fraproperties.newbounds = newbounds;

end