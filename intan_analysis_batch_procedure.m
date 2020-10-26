%% Intan Analysis Batch Procedure
% multiunit analysis for DMR, FRA, FRA10 and DMRREP
% Process recording data for each session/day
% file format: one .dat file for each channel (int16)
% first, .dat files are filtered and stored as .raw files (int16)
% then the voltage trace is thresholded to get mu activity, 
% which is then used to get mu response properties

% restoredefaultpath;
% addpath('E:\Congcong\Documents\MATLAB\MultiUnit_IntanAnalysis')
% addpath(genpath('E:\Congcong\Documents\MATLAB\support'))

probtype = {'H31x64', 'H22x32'};
%probtype = {'H31x64'};
stimuli = {'dmr', 'dmrrep', 'fra10', 'fra'};
%stimuli = {'dmr'};
flag_plot = 1;
flag_saveplot = 0;

sessions = dir('E:\Congcong\Documents\emsemble_thalamus\*CH');
%% data processing for each session
% STRF =  STA*#spike/pp/T; plot rfsig
% dmrrep: 2ms bin
% fra: smooth[2 2]
for session = 1:length(sessions)
    datafolder = ['E:\Congcong\Documents\emsemble_thalamus\', sessions(session).name];
    figpath = ['E:\Congcong\Documents\emsemble_thalamus\figure\multiunit\', sessions(session).name];
    if ~exist(figpath, 'dir')
        mkdir(figpath)
    end

    for ii = 1:length(stimuli)
        rawfolders = dir(fullfile(datafolder,sprintf('*_%s_*', stimuli{ii})));
        for jj = 1:length(rawfolders)
            cd(fullfile(datafolder, rawfolders(jj).name))
%             delete *.raw
%             delete *.mat
            badtrigfiles = intan_analysis_procedure(stimuli(ii), probtype, flag_plot, flag_saveplot, savepath);%
        end
    end
end





