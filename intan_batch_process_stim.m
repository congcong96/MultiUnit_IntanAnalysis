function outfile = intan_batch_process_stim(stimfolder, stimuli)
% intan_batch_process_acute_rat_tm Calculate FRAs, TM responses, and STRFs from multiunit data.
%
%    intan_batch_process_acute_rat_tm(stimfolder, stimuli) goes through all the directories
%    in an experimental folder and processes the data. This function should be
%    run inside a folder that contains individual folders for each penetration.
%    Inside the penetration folders are *-thresh.mat and *-trig.mat files, which
%    will be used to process the data.
%
%    stimfolder : absolute path to folder holding all the stimuli for the
%    experiments.
%       Example: stimfolder = 'F:\2013-9-18-trained-rat-tm20-acute\stimuli_2013-9-18';
%
%    stimuli : cell array holding strings of stimulus types that were presented
%
%    For experiments 2013-9-18, 2013-10-17, 2013-11-25:
%        stimuli =  {'fra', 'tm1',  'rn11',  'rn41',  'tm2',   'rn12',  'rn42'};
%
%    For experiment 2013-12-19:
%        stimuli =  {'fra1', 'fra2',  'tm', 'rn1',  'rn4'};
%
%    Note: you don't have to run all the stimuli at once. You can run
%    individual stimuli. In this case, you might set stimuli = {'fra1'},
%    for example.
%
%    Craig Atencio
%    12/5/13


%------------------------------------------------------
% Need function to convert intan channels to neuronexus
% channels to appropriate depth
%------------------------------------------------------

% added dmr, dmrrep, fra10
% -- updated by Congcong, 11/08/2019 

if ( nargin ~= 2 )
    error('You need two input arguments.');
end % (if)


if ( nargin == 2 )
    if ( isempty(stimfolder) )
        error('You must enter stimfolder input argument.');
    end
    
    if ( isempty(stimuli) )
        error('You must enter stimuli input argument.');
    end
end


% % Files that hold the FRA and TM stimulus data
% tcparamsfile = fullfile(stimfolder, ...
%    'acute_rat_tms_freq_resp_area_center6khz_range5oct_nreps1_fs96000hz_params.mat');
% tmparamsfile = fullfile(stimfolder, ...
%    'acute_rat_temporal_modulation_stim_40db_40TM_0SM_fs96000Hz_params.mat');
%
% % Names of the ripple envelope files
% rn1_sprfile = fullfile(stimfolder, 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr');
% rn4_sprfile = fullfile(stimfolder, 'rn4-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr');
%
% % Load and save the FRA parameters
% load(tcparamsfile, 'frequency', 'attenuation');
%
% % FRA function assumes a struct is passed to it. Put variables in a struct.
% tcparams.atten = abs(attenuation);
% tcparams.freq = frequency;
%
% % Load the TM parameters
% load(tmparamsfile, 'tmparams');



for j = 1:length(stimuli)
    
    stimtype = stimuli{j};
    
    threshfile = dir( sprintf('*-site*-*um-*db-%s-*-fs*-thresh.mat', stimtype) );
    trigfile = dir( sprintf('*-site*-*um-*db-%s-*-fs*-ADC-00.mat', stimtype) );
    
    
    if ( ~isempty(trigfile) && ~isempty(threshfile) && length(trigfile) == 1 )
        load(threshfile.name, 'thresh');
        load(trigfile.name, 'trigger');
        index = strfind(threshfile.name, '.mat');
        base_name = threshfile.name(1:index-1);
        
        
        if ( strcmp(stimtype,'fra') ) % tuning curve
            outfile = sprintf('%s-fra.mat', base_name);
            d = dir(outfile);
            
            if 1%( isempty(d) ) % We haven't previously processed the data ...
                
                tcparamsfile = fullfile(stimfolder, ...
                    'tms_freq_resp_area_center4khz_range6oct_nreps1_fs96000hz-params.mat');
                
                % Load and save the FRA parameters
                load(tcparamsfile, 'frequency', 'attenuation');
                
                % FRA function assumes a struct is passed to it. Put variables in a struct.
                tcparams.atten = abs(attenuation);
                tcparams.freq = frequency;
                
                tc = calculate_tuning_curve(thresh,trigger,tcparams);
                save(outfile, 'thresh', 'trigger', 'tcparams', 'tc');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif ( strcmp(stimtype,'fra1') ) % tuning curve
            outfile = sprintf('%s-fra1.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                
                tcparamsfile = fullfile(stimfolder, ...
                    'acute_rat_tms_freq_resp_area_center2khz_range4oct_nreps1_fs96000hz_params.mat');
                
                % Load and save the FRA parameters
                load(tcparamsfile, 'frequency', 'attenuation');
                
                % FRA function assumes a struct is passed to it. Put variables in a struct.
                tcparams.atten = abs(attenuation);
                tcparams.freq = frequency;
                
                tc = calculate_tuning_curve(thresh,trigger,tcparams);
                save(outfile, 'thresh', 'trigger', 'tcparams', 'tc');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
            
        elseif ( strcmp(stimtype,'fra2') ) % tuning curve
            outfile = sprintf('%s-fra2.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                tcparamsfile = fullfile(stimfolder, ...
                    'acute_rat_tms_freq_resp_area_center10khz_range4oct_nreps1_fs96000hz_params.mat');
                
                % Load and save the FRA parameters
                load(tcparamsfile, 'frequency', 'attenuation');
                
                % FRA function assumes a struct is passed to it. Put variables in a struct.
                tcparams.atten = abs(attenuation);
                tcparams.freq = frequency;
                outfile = sprintf('%s-fra2.mat', base_name);
                tc = calculate_tuning_curve(thresh,trigger,tcparams);
                save(outfile, 'thresh', 'trigger', 'tcparams', 'tc');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif ( strcmp(stimtype,'fra10') )
            outfile = sprintf('%s-fra10.mat', base_name);
            d = dir(outfile);
            if ( isempty(d) ) % We haven't previously processed the data ...
                tcparamsfile = fullfile(stimfolder, ...
                    'freq_resp_area_stimulus_flo500Hz_fhi32000Hz_nfreq21_natten8_nreps10_fs96000_param.mat');
                
                % Load and save the FRA parameters
                load(tcparamsfile, 'frequency', 'attenuation');
                
                % FRA function assumes a struct is passed to it. Put variables in a struct.
                tcparams.atten = abs(attenuation);
                tcparams.freq = frequency;
                outfile = sprintf('%s-fra10.mat', base_name);
                tc = calculate_tuning_curve(thresh,trigger,tcparams);
                save(outfile, 'thresh', 'trigger', 'tcparams', 'tc');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif contains(stimtype,'tm') % temporal modulation
            outfile = sprintf('%s-tm.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                tmparamsfile = fullfile(stimfolder, ...
                    'acute_rat_temporal_modulation_stim_40db_40TM_0SM_fs96000Hz_params.mat');
                load(tmparamsfile, 'tmparams');
                [tmresp] = calculate_temporal_modulation_response(thresh, trigger, tmparams);
                save(outfile, 'thresh', 'trigger', 'tmparams', 'tmresp');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif contains(stimtype,'rn1') % dynmaic moving ripple
            outfile = sprintf('%s-strf.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                rn1_sprfile = fullfile(stimfolder, 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr');
                strf = calculate_strf(thresh, trigger, 0, rn1_sprfile);
                save(outfile, 'thresh', 'trigger', 'strf');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif contains(stimtype,'rn4') % 4 DMRs added together
            outfile = sprintf('%s-strf.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                rn4_sprfile = fullfile(stimfolder, 'rn4-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr');
                strf = calculate_strf(thresh, trigger, 0, rn4_sprfile);
                save(outfile, 'thresh', 'trigger', 'strf');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif contains(stimtype,'rn8') % 8 DMRs added together
            outfile = sprintf('%s-strf.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                rn8_sprfile = fullfile(stimfolder, 'rn8-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr');
                strf = calculate_strf(thresh, trigger, 0, rn8_sprfile);
                save(outfile, 'thresh', 'trigger', 'strf');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif contains(stimtype,'rn16') % 16 DMRs added together
            outfile = sprintf('%s-strf.mat', base_name);
            d = dir(outfile);
            
            if ( isempty(d) ) % We haven't previously processed the data ...
                rn16_sprfile = fullfile(stimfolder, 'rn16-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr');
                strf = calculate_strf(thresh, trigger, 0, rn16_sprfile);
                save(outfile, 'thresh', 'trigger', 'strf');
                fprintf('Saved data in %s\n\n', outfile);
            else
                fprintf('Data exists in %s\n\n', outfile);
            end % if ( length(d)==0 )
            
        elseif contains(stimtype,'dmr') % dynmaic moving ripple
            
            if contains(stimtype,'rep')
                outfile = sprintf('%s-raster.mat', base_name);
                d = dir(outfile);
                
                if  isempty(d) % We haven't previously processed the data ...
                    for i = 1:length(thresh)
                        stimlen = 10000;
                        binsize = 2;
                        [raster_mat, spktimes_trial] = calc_dmrrep_raster(thresh(i), trigger, stimlen, binsize);
                        raster(i).exp = thresh(i).exp;
                        raster(i).chan = thresh(i).chan;
                        raster(i).probe = thresh(i).probe;
                        raster(i).stim = thresh(i).stim;
                        raster(i).atten = thresh(i).atten;
                        raster(i).raster_mat = raster_mat;
                        raster(i).spktimes_trial = spktimes_trial;
                        raster(i).fs = thresh(i).fs;
                        raster(i).stimlen = stimlen;
                        raster(i).binsize = binsize;
                    end
                    save(outfile, 'thresh', 'trigger', 'raster');
                    fprintf('Saved data in %s\n\n', outfile);
                else
                    fprintf('Data exists in %s\n\n', outfile);
                end
                
            else
                outfile = sprintf('%s-strf.mat', base_name);
                d = dir(outfile);
                
                if  isempty(d) % We haven't previously processed the data ...
                    dmr_sprfile = fullfile(stimfolder, 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506.spr');
                    %dmr_sprfile = fullfile(stimfolder, 'dmr-50flo-40000fhi-4SM-150TM-40db-96khz-96DF-20min.spr');
                    strf = calculate_strf(thresh, trigger, 0, dmr_sprfile);
                    save(outfile, 'thresh', 'trigger', 'strf');
                    fprintf('Saved data in %s\n\n', outfile);
                else
                    fprintf('Data exists in %s\n\n', outfile);
                end
            end
        end
    end % (if)
    
end % (for j)

return;







