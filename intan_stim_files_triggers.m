function badtrigfiles = intan_stim_files_triggers(stimuli, multprobeopt, udbandwidth)
%automatically_save_triggers Process Individual stim/trig files from long .ncs files
%
% automatically_save_triggers(stimtype) finds individual stimulus and
% trigger files withing larger .ncs files, and writes out each file to
% a separate smaller file. The function also finds the triggers within each
% smaller trigger file, and saves these to a .mat file.
%
% Note: The function must be inside the folder where the .Ncs files are located.
%
% stimtype describes the type of data. It is either 'ripples' or 'natural'.
% For 'ripples', only ripple stimuli were presented. For 'natural', a
% sequence of natural sounds was also presented. If stimtype is not
% included as an input argument, the default is 'ripples'.
%
% This function also makes the following function calls:
%
% ncs2raw_dir;
% [trigger] = findtrig_longfile(file,-0.6);
% extract_stim_trig_from_long_files(trigger, stimuli);
% [trigger, outfile] = findtrig(file, -0.6);
%
% Craig Atencio
% 2012-5-15
% 2012-7-11
% 2014-8-1
%automatically_save_triggers(stimtype);

badtrigfiles = {};

if strcmp(stimuli{1}, 'auto')
    stimuli = regexp(pwd, '(?<=(db_))\S+(?=(_\w+-\w+))','match','once');
end

trigfile = dir('*ADC-00.raw');

if ~isempty(trigfile)
    fprintf('\n%s already processed! Skipping...\n', pwd)
    pause(1)
    return
end

[amplifier_channels, board_adc_channels, frequency_parameters] = read_intan_rhd_file;
fs = frequency_parameters.board_adc_sample_rate;

if exist('udbandwidth', 'var')
    [trigfile, stim, outfile] = intan2raw(amplifier_channels, frequency_parameters, board_adc_channels, multprobeopt, udbandwidth);
else
    [trigfile, stim, outfile] = intan2raw(amplifier_channels, frequency_parameters, board_adc_channels, multprobeopt);
end


% Find the first trigger in the long trigger file
thresh = 0.6;
[trigger] = intan_findtrig_longfile(trigfile, thresh);
%trigger = trigger(1801:end);
if multprobeopt == 1 && strcmp(stimuli{:},'dmr')
    save([trigfile(1:end-3) 'mat'], 'trigger')
    return %modifield for simultaneous recording
end

% Get the individual stimulus and trigger .raw files
% for spontaneous activities, just rename the file to include min
if strcmp(stimuli{:}, 'spon')
    files = dir(sprintf('*-%s-*-A-*.raw', 'spon')); 
    for n = 1:length(files)
        infile = files(n).name;
        part1exp = '^\S+(?=(-spon))';
        part2exp = '(?<=(spon-))\S+$';
        part1 = regexp(infile, part1exp, 'match', 'once');
        part2 = regexp(infile, part2exp, 'match', 'once');
        outfile = sprintf('%s-%s-%s-%s', part1, 'spon', '15min', part2);
        movefile(infile,outfile);
    end
    return
end
[~] = intan_extract_stim_trig_from_long_files(trigger, fs, stim, stimuli);

% if strfind(stimuli{:}, 'rep')
%     trigger = first_triggers;
%     d = dir(sprintf('*-%s-*-ADC-*.raw', stimuli{1}));
%     if isempty(d)
%         d = dir('*-ADC-*.raw');
%     end
%     file = d.name;
%     outfile = [file(1:end-10) 'ADC-00.mat'];
%     if length(unique(diff(trigger)))>2
%         warning('Something wrong with %s trigger',outfile)
%         badtrigfiles = [badtrigfiles, outfile];
%     end
%     
%     save(outfile,'trigger');
%     clear('file', 'trigger', 'outfile');
%     return
% end
% delete original raw files if not required (e.g. when only one stimulus is
% present per folder)
if length(stimuli) == 1
    cellfun(@delete, outfile);
    delete(trigfile);
end


% Get the triggers from the trigger *.raw files
for i = 1:length(stimuli)
   d = dir(sprintf('*-%s-*-ADC-*.raw', stimuli{i})); 
   file = d.name;
   [trigger, outfile] = findtrig(file, thresh);
   trigger = trigger(trigger>20000-1);
   if length(unique(diff(trigger)))>2
       warning('Something wrong with %s trigger',outfile)
       badtrigfiles = [badtrigfiles, outfile];
   end

   save(outfile,'trigger');
   clear('file', 'trigger', 'outfile');
end % (for i)



return;



