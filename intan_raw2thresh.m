%THRESHOLD2MAT get threshold crossing for .raw files in given directory
%   [THRESH,OUTFILE] = THRESHOLD2MAT(DPATH,NSTD,STIMTYPE,FS) similar in
%   function to spikesort2mat, this function returns threshold crossing
%   times instead of spike times.  NSTD determines the threshold level in
%   number of recording standard deviations.
%
%   Written by Jonathan Shih 2010-11-10
function [thresh,outfile] = intan_raw2thresh(dpath, nstd, stimtype)

%Checking input parameters
if ( nargin ~= 3 )
    error('You need 3 input args.');
end


%Creating THRESH structs
str = struct(...
    'file',        [], ...
    'exp',         [], ...
    'site',        [], ...
    'amplifier',   [], ...
    'chan',        [], ...
    'depth',       [], ...
    'stim',        [], ...
    'atten',       [], ...
    'spiketimes',  [], ...
    'fs',          []);


%Removing backslash at end of directory path if necessary
if ( length(dpath) > 0 )
   if(dpath(end) == '\')
       dpath = dpath(1:(end-1));
   end
end


% Example file: 140910-site1-1606um-30db-rn1-fs20000-A-000.raw
dirfiles = gfn(fullfile(dpath, sprintf('*-%s-*min*.raw', stimtype)), 1);

if isempty(dirfiles)
    dirfiles = gfn(fullfile(dpath, sprintf('*-%s-*.raw', stimtype)), 1);
end

% Don't want the trigger channels, which are 'ADC' files
trigidx = ~contains(dirfiles, 'ADC');
datafiles = dirfiles(trigidx);

thresh = [];
%Going through each channel to get the threshold crossing times
ii = length(datafiles);
rawfile = datafiles{ii};
[~,b,c] = fileparts(rawfile);
basefile = [b c];
basename = regexp(basefile, '\S+(?=(-\w{1}-\d{3}.raw))','match','once');
outfile = fullfile(dpath, [basename '-thresh.mat']);
d = dir(outfile);
if ~isempty(d)
    thresh = [];
    return
end
for ii = 1:length(datafiles)
    %Current file being processed
    rawfile = datafiles{ii};
    [~,b,c] = fileparts(rawfile);
    basefile = [b c];
    
    details = regexp(basefile, ['(?<exp>\d{6}_\d{6})-site(?<site>\d{1,2})-'...
        '(?<depth>\d{3,})um-(?<atten>\w+db)-' sprintf('(?<stim>%s)-',stimtype)...
        '(?<stimlen>\d{1,2})min-(?<probe>\S+)-fs'...
        '(?<fs>\d{3,6})-(?<amplifier>[ABCD])-(?<chan>\d{3}(?=(.raw)))'],...
        'names');
    
    if isempty(details)
        
        details = regexp(basefile, ['(?<exp>\d{6}_\d{6})-site(?<site>\d{1,2})-'...
            '(?<depth>\d{3,})um-(?<atten>\w+db)-' sprintf('(?<stim>%s)-',stimtype)...
            '(?<probe>\S+)-fs'...
            '(?<fs>\d{3,6})-(?<amplifier>[ABCD])-(?<chan>\d{3}(?=(.raw)))'],...
            'names');
        
    end
    
    if isempty(details)
        
        details = regexp(basefile, ['(?<exp>\d{6}_\d{6})-site(?<site>\d{1,2})-'...
            '(?<depth>\d{3,})um-(?<atten>\w+db)-' sprintf('(?<stim>%s)-',stimtype)...
            'fs(?<fs>\d{3,6})-(?<amplifier>[ABCD])-(?<chan>\d{3}(?=(.raw)))'],...
            'names');
        
    end

   
    fs = str2double(details.fs);
    fidin = fopen(rawfile, 'r');
    nsamples = 60 * fs;
    [samples, ~] = fread(fidin, nsamples, 'int16');
    fclose(fidin);
    
    %Getting STD of .raw file
    stdsample = std(samples);
    
    %Initializing struct
    s = str;
    s.file = basefile;
    s.exp = details.exp;
    s.site = str2double(details.site);
    s.amplifier = details.amplifier;
    s.depth = str2double(details.depth);
    s.chan = str2double(details.chan); % assign channel number
    s.stim = stimtype;
    if isfield(details, 'stimlen')
        s.stimlen = str2double(details.stimlen);
    end
    s.atten = details.atten;
    if isfield(details, 'probe')
        s.probe = details.probe;
    else
        s.probe = 'a1x32-poly3';
    end
    s.spiketimes = get_threshold_crossings(rawfile,nstd*stdsample,fs);
    s.fs = fs; % save sampling frequency
    

    %Updating user
    fprintf('Finished %s\n', basefile);
    
    thresh = [thresh s];
    clear('s');
end




% outfile = sprintf('%s-site%d-%dum-%ddb-%s-fs%d-thresh.mat', ...
%    exp, site, depth, atten, stimtype, fs);

return;




%Helper Function: get_threshold_crossings
%----------------------------------------
%Gets the threshold crossings from a .raw file and converts them into msec
function crossings = get_threshold_crossings(rawpath,thresh,fs)

%Opening .raw file
fid = fopen(rawpath,'r');

%Number of data values to read in each iteration
blockSize = 10*1024;

init_value = 0;
nblocks = 0;
crossings = [];
while(~feof(fid))
    %Getting current block of data
    blockData = fread(fid,blockSize,'int16');
    blockData = -blockData; % invert waveform because large spike peaks are downward in Intan
                            % software; they are upward in neuralynx
                            % software
    
    %Appending last value of previous block to current block
    blockData = [init_value; blockData];
    
    %Finding threshold crossings in the current block
    ind_supra = find(blockData(2:end) > thresh);
    ind_cross = ind_supra(find(blockData(ind_supra) <= thresh));
    crosstimes = (nblocks*blockSize + ind_cross - 1)/fs*1000; %converting samples into msec
    crossings = [crossings; crosstimes];
    
    %Updating value to append to beginning of next block
    init_value = blockData(end);
    
    %Updating number of blocks
    nblocks = nblocks + 1;
    
%     plot(blockData), hold on;
%     plot(xlim,thresh*[1 1],'k--');
%     scatter(ind_cross,blockData(ind_cross+1),'rx'), hold off;
%     pause;
end

%Closing .raw file
fclose(fid);

return;



