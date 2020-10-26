function raster = calculate_raster(spk,trigger,params)

% calculate_raster(spk,trigger,params)
% modified calculate tuning curve for raster plot
% OUTPUT
% raster has spiketimes for each presentation sorted by atten and freq
% raster.atten = vector of stimulus attenuations for raster.rastermat
% raster.freq  = vector of stimulus frequencies for raster.rastermat
% raster.rastermat = matrix of post-stimulus spike counts for each stimulus presentation
%
% c.f.
% CALCULATE_TUNING_CURVE calculate tuning curve
%   TC = CALCULATE_TUNING_CURVE(SPK,TRIGGER,PARAMS) calculates the tuning
%   curves given the spike times stored in SPK and the trigger times stored
%   in TRIGGER.  TC has the following struct fields:
%
%       TC.atten = vector of stimulus attenuations
%       TC.freq  = vector of stimulus frequencies
%       TC.tcmat = matrix of post-stimulus spike counts for each stimulus
%                  condition
%
%   TC = CALCULATE_TUNING_CURVE(SPK,TRIGGER,PARAMS,TLIM) allows the user to
%   specify the post-stimulus time window (in msec) in which spikes will be
%   included in the spike counts. (default: TLIM = [0 50])
%
%   Written by Jonathan Shih 5-4-2007
%   check if 'atten' exists as a field of spk, or named as 'attenuation'
%   raster(kk).rastermat(rasterInd,1) = {[raster(kk).rastermat{rasterInd,1} spktimes]};
%   to accomadate to repeats
%   updated by Congcong Hu, 12-20-2019


%Checking input arguments
if(nargin < 3)
    error('Need at least 3 input arguments.');
end
if isfield(params, 'atten')
    if(length(params.atten) ~= length(trigger))
        error('Number of triggers and number of stimulus conditions much match');
    end
else
    if(length(params.attenuation) ~= length(trigger))
        error('Number of triggers and number of stimulus conditions much match');
    end
end

%Converting trigger times into msec (from samples)
if ~isfield(spk,'fs')
    triggerTime = 1000*trigger/20000;
else
    triggerTime = 1000*trigger/spk(1).fs;
end
%Creating stimulus condition vectors
if isfield(params, 'atten')
    attenVec = sort(unique(params.atten));
    freqVec = sort(unique(params.freq));
else
    attenVec = sort(unique(params.attenuation));
    freqVec = sort(unique(params.frequency));
end
attenVec2 = repmat(attenVec,[length(freqVec),1]);
freqVec2 = repmat(freqVec,[1,length(attenVec)]);
freqVec2 = freqVec2';
freqVec2 = freqVec2(:);

%Calculate the difference of triggers
Dtrigger = diff(trigger);
MinDtrigger = min(Dtrigger); % c.f. unique(Dtrigger) => 6000 and 6001
if ~isfield(spk,'fs')
    tlim = [0 1000*MinDtrigger/20000]; % in ms
else
    tlim = [0 1000*MinDtrigger/spk(1).fs]; % in ms
end
%Initializing raster struct
for ii = 1:length(spk)
    raster(ii).Dtrigger = tlim(2);
    raster(ii).atten = attenVec2;
    raster(ii).freq = freqVec2;
    raster(ii).rastermat = cell(length(attenVec2),1);
end

%Finding spikes
for jj = 1:length(trigger)
    %Getting parameters for current trial
    if isfield(params, 'atten')
        atten = params.atten(jj);
        freq = params.freq(jj);
        stimTime = triggerTime(jj);
    else
        atten = params.attenuation(jj);
        freq = params.frequency(jj);
        stimTime = triggerTime(jj);
    end
    
    %Finding indices in raster matrix for trial stimulus conditions
    rasterInd = intersect(find(attenVec2 == atten),find(freqVec2 == freq));
    
    for kk = 1:length(spk)
        spktimesInd =  (spk(kk).spiketimes >= (stimTime + tlim(1))) & (spk(kk).spiketimes <= (stimTime + tlim(2))) ;
        spktimes = spk(kk).spiketimes(spktimesInd) - stimTime;
        spktimes = spktimes(:);
        raster(kk).rastermat(rasterInd,1) = {[raster(kk).rastermat{rasterInd,1}; spktimes]};%mat2cell(spktimes,length(spktimes),1)
    end
end