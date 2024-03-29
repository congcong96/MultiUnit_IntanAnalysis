function [sig, RI, z, strfcorr, nstrfcorr] = calc_strf_RI(stimulus, spk, trigger, niter)
% Based on Escabi et al., 2014 (J. Neurophysiol)
% 
% Determine if the STRF is stable over time
% The RI of the STRF is compared to the null distribution of RIs of nSTRFs
%
% Inputs:
%   stimulus: file path of .spr file
%   spk: spk structure of individual neuron.
%   trigger: trigger of the recording
%   niter: number of iteration to calculate reliability index distribution, 
%             default = 1000
%   
%
% Outputs:
%   sig: significance of the reliability index of the unit
%   RI: reliability index of the unit
%   strfcorr: (RI_iter * 1) distribution of correlation of the neuron
%   nstrfcorr: (RI_iter * 1) null distribution of correlation
% null correlation distribution is generated by circular shift spike train once and sample spikes
% for RI_niter times
%
% Written by Congcong, 04/01/2020.
    

% RELIABILITY INDEX OF STRF

if nargin == 3
    niter = 100;
end
stimstart = trigger(1)/20;
stimend = (trigger(end)+trigger(2)-trigger(1))/20;
dur = stimend-stimstart;%in ms
spiketimes = spk.spiketimes;
spiketimes(spiketimes<stimstart) = [];
spiketimes(spiketimes>stimend) = [];
n = floor(length(spiketimes)/2);
strfcorr = zeros(niter,1);
for ii = 1:niter
     if mod(ii, 100) == 0
        fprintf('interation %d/%d\n', ii, niter)
    end
    newtrain = spiketimes(randperm(length(spiketimes)));
    spk.spiketimes = newtrain(1:n);
    strf1 = calculate_strf_v2(spk, trigger, stimulus, 0.2, 0.05);
    spk.spiketimes = newtrain(n+1:end);
    strf2 = calculate_strf_v2(spk, trigger, stimulus, 0.2, 0.05);
    strfcorr(ii) = corr(strf1.rfcontra(:), strf2.rfcontra(:));
    %strfcorr(ii,2) = corr(strf1.rfipsi(:), strf2.rfipsi(:));
    spk.spiketimes = spiketimes;
end

spiketimes = spiketimes + rand * dur;
idx = find(spiketimes > stimend, 1);
if ~isempty(idx)
    if idx > 1
        spiketimes= [spiketimes(idx:end)-stimend; spiketimes(1:idx-1)] ;
    else
        spiketimes= spiketimes - stimend;
    end
end
nstrfcorr = zeros(niter,1);
for ii = 1:niter
    if mod(ii, 100) == 0
        fprintf('interation %d/%d\n', ii, niter)
    end
    newtrain = spiketimes(randperm(length(spiketimes)));
    spk.spiketimes = newtrain(1:n);
    strf1 = calculate_strf_v2(spk, trigger, stimulus, 0.2, 0.05);
    spk.spiketimes = newtrain(n+1:end);
    strf2 = calculate_strf_v2(spk, trigger, stimulus, 0.2, 0.05);
    nstrfcorr(ii) = corr(strf1.rfcontra(:), strf2.rfcontra(:));
    %nstrfcorr(ii,2) = corr(strf1.rfipsi(:), strf2.rfipsi(:));
    spk.spiketimes = spiketimes;
end
RI = mean(strfcorr) - mean(nstrfcorr);
strfsem = std(strfcorr)^2/niter;
nstrfsem = std(nstrfcorr)^2/niter;
z = RI/sqrt(strfsem+nstrfsem);
sig =  z > 2.33;

