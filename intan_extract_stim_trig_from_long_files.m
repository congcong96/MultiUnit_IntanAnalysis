function first_triggers = intan_extract_stim_trig_from_long_files(trigger, fs, folderstim, stimuli)
% extract_stim_trig_from_long_files Get individual files from within a long file
%
%   extract_stim_trig_from_long_files(trigger, stimuli)
%
%   This function is meant to accompany the data for acute rat temporal
%   modulation experiments.
%
%   This function uses a template of first triggers, i.e. what is the 
%   spacing between the first trigger for each stimulus. Using this 
%   template, the files are extracted.
%
%   trig1 : first trigger in long trigger file. Supplied by calling
%   function.
%
%   stimtype : type of stimuli. May be either 'allstim' or may be a list of
%   stimuli. The total list of stimuli is
%
%   totalstimuli =  {'fra', 'tm1',  'rn11',  'rn41',  'tm2',   'rn12',  'rn42'};
%
%   If the files contain less than this, then the index to the stimuli
%   should be specified. For example, if only 'fra'  and 'tm1' stimuli were
%   presented, then stimtype = 1:2. If 'rn41' and 'tm2' stimuli were
%   presented, then stimtype = 4:5.
%
%
% caa 9/20/13
%
% extract_stim_trig_from_long_files(template, stimuli);
% return first_triggers
% if repeat stimuli, return without truncating .raw files
% -- updated by Congcong, 11/08/2019


% Parse the input arguments
narginchk(4,4)

if isempty(stimuli)    
    stimuli = {folderstim};
end


% Now find where the stimuli and triggers are within the larger files:
d = diff(trigger);
index = find(d > 9 * fs); % Search for deadtimes greater than 9 seconds
                          % There should be at least 9 seconds b/w stimuli 
first_triggers = [trigger(1) trigger(index+1)]; % #first triggers = number of stimuli
% first_triggers = first_triggers(1:length(stimuli));

% Make sure there wasn't a bug %temp
if ( length(first_triggers) ~= length(stimuli) ) && (length(first_triggers) ~= length(trigger))
   error('Number of first triggers and number of stimuli do not match.');
end

if length(first_triggers) == length(trigger) % for dmrrep, triger intervals are larger than 9s
    first_triggers = first_triggers(1);
end

% The number of samples to save before the start of the stimulus
delta = 1*fs; % Get 1 sec of deadtime before each individual stimulus starts

% How many samples to process at a time
nsamples = 1024 * 256;


% Get the raw files within the present directory; These include the
% amplifier channel files and the trigger file. These files hold the data
% for all stimuli. Each individual file will be decomposed into multiple
% smaller files, with each one holding data for one of the stimuli.
files = dir( sprintf('*-site*-*um-*db-%s-*.raw',folderstim) ); % length(files) == 17: 16 channels, 1 trigger

% Now go through each long file, and extract the shorter stimulus files
% within each longer file; this will work for the trigger files too, since
% it goes through every *.raw file
for n = 1:length(files)

   % File that we will process in this interation
   infile = files(n).name;   

   % Get the components of the trigger file; we'll need these when we 
   % write out the smaller files, since each smaller file will have to be
   % correctly named.
   
   part1exp = ['^\S+(?=(-' folderstim '))'];
   part2exp = ['(?<=(' folderstim '-))\S+$'];
   part1 = regexp(infile, part1exp, 'match', 'once');
   part2 = regexp(infile, part2exp, 'match', 'once');

   
   for i = 1:length(stimuli) % go through each smaller file within the larger file
                                    % length(first_triggers) == number of smaller files
      start = first_triggers(i) - delta; % the starting sample number

      % The ending sample number
      if ( i < length(first_triggers) )
         finish = first_triggers(i+1) - delta;
      else % it's the last stimulus
         finish = trigger(end) + fs ;
      end
      
      stimlen = round((finish-start)/fs/60); %length of stimulus in minutes
      stimlen = [num2str(stimlen, '%d') 'min'];

      % Make the output file name
      outfile = sprintf('%s-%s-%s-%s', part1, stimuli{i}, stimlen, part2);
      
      if exist(outfile,'file') == 2 %isempty(d) )
         
          fprintf('%s already exists\n\n', outfile);
          
      else
          fprintf('\n%s\n%s\n[start = %.0f, finish = %.0f]\n', infile, outfile, start, finish);
          
          % Extract the individual files within the larger file
          fidin = fopen(infile,'r');
          fidout = fopen(outfile, 'w');

          ii = 1;

          while ( ~feof(fidin) )

             [data, count] = fread(fidin, nsamples,'int16');
             index = ((ii-1)*nsamples+1):(ii*nsamples);

             % Get the indices to the data between the triggers
             ind_good = find( index >= start & index <= finish );

             if ( count < length(ind_good) )
                ind_good = ind_good(1:count);
             end

             fwrite(fidout, data(ind_good), 'int16');

             ii = ii + 1;

          end

          fclose(fidin);
      
      
      end % (if isempty(d))

   end % (for i)

   fclose('all');

end % (for n)

fclose('all');

return;


