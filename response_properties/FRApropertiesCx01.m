function [fraproperties]=FRApropertiesCx01(bounds,freqs,spikes,levels)
%% FRAproperties05
% this script is for calculating threshold, BF, CF, Q10, Q40, ratio of tuning curve (Tratio) and peaks of tuning curve (NumPeak)
% woking with following the script drawFRAbounds03
% original script is FRApropertiesCx01

%% extract BF (best frequency) - a frequency with the maximum response
% nfreqs=size(unique(freqs),1);
% freqspace=linspace(min(log(freqs)),max(log(freqs)),nfreqs*2);
% maxspikes=max(max(spikes));
% [f1 f2]=find(spikes==maxspikes);
% BFvals=mean(freqspace(f2*2));
% BF=exp(BFvals)/1000;

%% extract BF2 (another dbest frequency) - a frequency with the maximum response across all levels
nfreqs=size(unique(freqs),1);
freqspace=linspace(min(log(freqs)),max(log(freqs)),nfreqs);
acrosslevels=sum(spikes,1);
[~, ff]=max(acrosslevels);
BFvals2=mean(freqspace(ff));
BF2=exp(BFvals2)/1000;

%%  if there are frequency gaps in the lowest intensity, take a part that has BF.
[I2]=find(bounds<=(length(unique(levels)))); % find the bound in FRA
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
if min(newbounds)>=length(unique(levels))
    CF = NaN;
    BW10freq = [];
    Q10=NaN;
    BW30freq = [];
    Q30 = NaN;
    BW40freq = [];
    Q40 = NaN;
else
    [~, cf] = min(newbounds);
    if ~isempty(diff(cf)) && ~any(diff(cf)~=1)
        cf = mean(cf); %if the tip is 'dull', take the average;
    end
    
    if length(cf)>1
        [~, idx] = min(abs(cf-bf)); %if multiple peak, take CF closest to BF
        cf = cf(idx);
    end
    CF = exp((freqspace(ceil(cf))+freqspace(floor(cf)))/2)/1000;
    
    %get the bandwidth for Q10
    if (levels(2) - levels(1)) == 10
        q=(min(newbounds)+1);
    elseif (levels(2) - levels(1)) == 5
        q=(min(newbounds)+2);
    end
    I = newbounds <= q;
    
    bw=floor([cf cf]);
    while bw(1) > 1 && I(bw(1)-1) %moveleft until hit edge of I
        bw(1) = bw(1) - 1;
    end
    
    while bw(2) < nfreqs && I(bw(2)+1) %moveleft until hit edge of I
        bw(2) = bw(2) + 1;
    end
    
    highfreq=exp(freqspace(bw(2)))/1000;
    lowfreq=exp(freqspace(bw(1)))/1000;
    
    BW10freq(2)=highfreq;
    BW10freq(1)=lowfreq;
    bandwidth= exp((diff(bw)+1)*(freqspace(2)-freqspace(1)))/1000;
    Q10=CF/bandwidth;
    
    %% Extract Q30 (repeate for Q30)
    % Q30 = CF / the width of the tuning curve 30dB above threshold
    
    if (levels(2) - levels(1)) == 10
        q=(min(newbounds)+3);
    elseif (levels(2) - levels(1)) == 5
        q=(min(newbounds)+6);
    end
    if q <= length(levels)
        [I]=find(newbounds<=q);
        bw=floor([cf cf]);
        while bw(1) > 1 && I(bw(1)-1) %moveleft until hit edge of I
            bw(1) = bw(1) - 1;
        end
        
        while bw(2) < nfreqs && I(bw(2)+1) %moveleft until hit edge of I
            bw(2) = bw(2) + 1;
        end
        
        highfreq=exp(freqspace(bw(2)))/1000;
        lowfreq=exp(freqspace(bw(1)))/1000;
        
        BW30freq(2)=highfreq;
        BW30freq(1)=lowfreq;
        bandwidth= exp((diff(bw)+1)*(freqspace(2)-freqspace(1)))/1000;
        Q30=CF/bandwidth;
        
        %% Extract Q40 (repeate for Q40)
        if (levels(2) - levels(1)) == 10
            q=(min(newbounds)+4);
        elseif (levels(2) - levels(1)) == 5
            q=(min(newbounds)+8);
        end
        if q <= length(levels)
            [I]=find(newbounds<=q);
            bw=floor([cf cf]);
            while bw(1) > 1 && I(bw(1)-1) %moveleft until hit edge of I
                bw(1) = bw(1) - 1;
            end
            
            while bw(2) < nfreqs && I(bw(2)+1) %moveleft until hit edge of I
                bw(2) = bw(2) + 1;
            end
            
            highfreq=exp(freqspace(bw(2)))/1000;
            lowfreq=exp(freqspace(bw(1)))/1000;
            
            BW40freq(2)=highfreq;
            BW40freq(1)=lowfreq;
            bandwidth= exp((diff(bw)+1)*(freqspace(2)-freqspace(1)))/1000;
            Q40=CF/bandwidth;            
            
        else
            Q40=NaN;
            BW40freq=[];
        end
        
        
    else
        Q30=NaN;
        BW30freq=[];
        Q40=NaN;
        BW40freq=[];
    end
end
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
fraproperties.newbounds = newbounds;