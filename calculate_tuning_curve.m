%CALCULATE_TUNING_CURVE calculate tuning curve
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

function tc = calculate_tuning_curve(spk,trigger,params,tlim)

%Checking input arguments
if nargin < 3
    error('Need at least 3 input arguments.');
end

if isfield(params, 'atten') && (length(params.atten) ~= length(trigger))
    error('Number of triggers and number of stimulus conditions much match');
elseif isfield(params, 'attenuation') && (length(params.attenuation) ~= length(trigger))
end

%Setting default value for post-stimulus time window
if(~exist('tlim','var'))
    tlim = [0 50];
end

%Converting trigger times into msec (from samples)
if ~isfield(spk, 'fs')
    trigger = 1000*trigger/20000;
else
    trigger = 1000*trigger/spk(1).fs;
end

%Creating stimulus condition vectors
if isfield(params, 'atten')
    attenVec = sort(unique(params.atten));
    freqVec = sort(unique(params.freq));
else
    attenVec = sort(unique(params.attenuation));
    freqVec = sort(unique(params.frequency));
end

if attenVec(1) >= 0
    attenVec = sort(-attenVec);
end

%Initializing tuning curve struct
for ii = 1:length(spk)
    tc(ii).atten = attenVec;
    tc(ii).freq = freqVec;
    tc(ii).tcmat = zeros(length(attenVec),length(freqVec));
    if isfield(spk, 'probe')
        tc(ii).probe = spk(ii).probe;
        tc(ii).amplifier = spk(ii).amplifier;
    end
    tc(ii).chan = spk(ii).chan;
    
end

%Finding spike counts
for jj = 1:length(trigger)
    %Getting parameters for current trial
    if isfield(params, 'atten')
        atten = params.atten(jj);
        freq = params.freq(jj);
    else
        atten = params.attenuation(jj);
        freq = params.frequency(jj);
    end
    stimTime = trigger(jj);

    %Finding indices in tuning curve matrix for trial stimulus conditions
    tcattenInd = find(abs(attenVec) == abs(atten));
    tcfreqInd = find(freqVec == freq);

    for kk = 1:length(spk)
        spkCount = length(find( (spk(kk).spiketimes >= (stimTime + tlim(1))) & (spk(kk).spiketimes <= (stimTime + tlim(2))) ));
        tc(kk).tcmat(tcattenInd,tcfreqInd) = tc(kk).tcmat(tcattenInd,tcfreqInd) + spkCount;
    end
end