function [raster, spktimes_trial] = calc_dmrrep_raster(nevstruct, trigger, stimlen, binsize)

nreps = length(trigger);
spktimes = nevstruct.spiketimes;
if ~isfield(nevstruct, 'fs')
    fs = 20000;
else
    fs = nevstruct.fs;
end

trigger_time = trigger/fs*1000;
dtrigger = trigger_time(2) - trigger_time(1);
if size(trigger_time, 1) == 1
    trigger_time = trigger_time';
end
trigger_time(:,2) = [trigger_time(2:end,1); trigger_time(end,1) + dtrigger];

spktimes_trial = cell(nreps,1);
for i = 1:length(spktimes_trial)
    spktimes_trial{i} = spktimes(spktimes >= trigger_time(i,1) & spktimes <= trigger_time(i,2)) - trigger_time(i,1);
end

edges = 0:binsize:stimlen;

raster = cell2mat(cellfun(@(x) histcounts(x, edges), spktimes_trial, 'UniformOutput', 0));
end