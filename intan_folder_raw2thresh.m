function intan_folder_raw2thresh(stimuli)
% intan_batch_raw2thresh Get multi-unit data from raw files in experiment folders
% 
%    intan_batch_raw2thresh is run outside of the data folders holding all the .raw
%    files. It expects that within each folder there are .raw files for each
%    type of stimulus. It then processes each .raw file, finds multi-unit
%    data, and saves the data. It then moves on to th next folder.
% 
%    The stimuli that this function looks for have specific labels in the file
%    name. 
%
%    For experiments 2013-9-18, 2013-10-17, 2013-11-25:
%        stimuli =  {'fra', 'tm1',  'rn11',  'rn41',  'tm2',   'rn12',  'rn42'};
% 
%    For experiment 2013-12-19:
%        stimuli =  {'fra1', 'fra2',  'tm', 'rn1',  'rn4'};
%
%    If you have different stimuli, then you can run batch_raw2thresh(stimuli),
%    where stimuli is a cell array of strings holding different stimulus types.
%
%    This function expects there to be 16 channels of data. If there are
%    more or less then the function will need to be modified.
% 
%    Craig Atencio 11/20/13


if ( nargin == 0 )
   error('You need to input the stimuli that were used.');
end

p = pwd;

dpath = '.';
nstd = 4;

   for j = 1:length(stimuli)

      stimtype = stimuli{j};

      stimfiles = dir( sprintf('*-%s-*-*-*.raw', stimtype) );


      if ( ~isempty(stimfiles) )

         fprintf('MU for %s, %s\n', p, stimtype);

        [thresh,outfile] = intan_raw2thresh(dpath, nstd, stimtype);
         d = dir(outfile);

         if ( isempty(d) )

            save(outfile, 'thresh');
            fprintf('Saved data in %s.\n\n', outfile);
         else
            fprintf('%s already exist.\n\n', outfile);
         end

      end % (if)

   end % (for j)

return;




