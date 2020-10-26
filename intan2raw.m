function [trig_outfile, stim, outfile] = intan2raw(amplifier_channels, frequency_parameters, board_adc_channels, multprobeopt,udbandwidth)
% intan2raw Convert Intan AD files to .raw files
% 
%    trig_outfile = intan2raw(amplifier_channels, frequency_parameters, board_adc_channels)
%    amplifier_channels : struct holding amplifier channel recording information
%    frequency_parameters : frequency specs for amplifier channel recording
%    board_adc_channels : struct holding external analog input information
% 
%    The input arguments are obtained using:
% 
%    [amplifier_channels, board_adc_channels, frequency_parameters] = ...
%       read_intan_rhd_file;
% 
%    The intan file that is read is info.rhd. intan2raw assumes that the
%    recordings were made in one-channel-per-file mode.
% 
%    intan2raw goes through each intan .data amplifier recording and converts
%    it to a signed 16 bit integer file. It then filters that file using 
%    the lowpass and highpass frequencies listed in frequency_parameters.
% 
%    The board_adc_channels struct holds the trigger channel input, and this
%    is also converted from an unsigned 16 bit integer file to a signed 16
%    bit integer file.
%
%    Updated by JS 6/28/17 to include probename in .raw file names


trig_native_channel_name = board_adc_channels.native_channel_name;
trigfile = sprintf('board-%s.dat', trig_native_channel_name);

fs = frequency_parameters.amplifier_sample_rate;
if exist('udbandwidth','var')
    lower_cutoff = udbandwidth(1);
    upper_cutoff = udbandwidth(2);
else
    lower_cutoff = frequency_parameters.desired_lower_bandwidth;
    upper_cutoff = frequency_parameters.desired_upper_bandwidth;
end

[~, folder] = fileparts(pwd);
if contains(folder, 'ECoG')
    details = regexp(folder, ['(?<site>^site\d+)_(?<depth>\S+um)_'...
        '(?<atten>\S+db)_(?<stim>\w{3,}?(?=[_E]))_(?<probe>[E]\S+)_'...
        '(?<exp>\d{6}_\d{6}$)'], 'names');
else
    details = regexp(folder, ['(?<site>^site\d+)_(?<depth>\S+um)_'...
        '(?<atten>\S+db)_(?<stim>\w{3,}?(?=[_Ha\d]))_(?<probe>[Ha]\d\S+)_'...
        '(?<exp>\d{6}_\d{6}$)'], 'names');
end

if isempty(details)
 details = regexp(folder, ['(?<site>^site\d+)_(?<depth>\S+um)_'...
     'hp(?<atten>\S+db)_*_(?<stim>\w{3,})_*_(?<exp>\d{6}_\d{6}$)'], 'names');
 details.probe = 'H31x64';
 details.stim = 'dmr';
end

if multprobeopt
    details.sep_depth = regexp(details.depth, '_', 'split');
    details.sep_probe = regexp(details.probe, '(?<=x\d{2,})_','split');
    uniqueamps = unique({amplifier_channels.port_prefix});
    assert(length(uniqueamps) == length(details.sep_probe))
    assert(length(uniqueamps) == length(details.sep_depth))
end

nsamples = 1024*256;

outfile = cell(length(amplifier_channels), 1);

for i = 1:length(amplifier_channels)

   native_channel_name = amplifier_channels(i).native_channel_name;
   amp_file = sprintf('amp-%s.dat', native_channel_name);
   
   if multprobeopt && strcmp(details.stim, 'dmr')
       mpidx = find(strcmp(uniqueamps, native_channel_name(1)));
       outfile{i} = sprintf('%s-%s-%s-%s-%s-30min-%s-fs%.0f-%s.raw', ...
           details.exp, details.site, details.sep_depth{mpidx}, details.atten,...% -30min-for simultaneous recording
           details.stim, details.sep_probe{mpidx}, fs, native_channel_name);
   elseif multprobeopt
       mpidx = find(strcmp(uniqueamps, native_channel_name(1)));
       outfile{i} = sprintf('%s-%s-%s-%s-%s-%s-fs%.0f-%s.raw', ...
           details.exp, details.site, details.sep_depth{mpidx}, details.atten,...
           details.stim, details.sep_probe{mpidx}, fs, native_channel_name);
   else
       outfile{i} = sprintf('%s-%s-%s-%s-%s-%s-fs%.0f-%s.raw', ...
           details.exp, details.site, details.depth, details.atten, details.stim,...
           details.probe, fs, native_channel_name);
   end
       
   fprintf('\n\n%s\n', outfile{i});
  
  
   d = dir(outfile{i});
   
   if ( isempty(d) )
       fprintf('%s\n', amp_file);


       %Reading and writing raw waveform data
       fpathin = amp_file;
       fpathout = outfile{i}; % one file change to use the addition

       fidin = fopen(fpathin,'r');
       fidout = fopen(fpathout,'w');

       while ~feof(fidin)
          data = fread(fidin,nsamples,'int16');
          if( ~isempty(data) )
             fwrite(fidout,data,'int16');
          end
       end % (while)

       fclose('all');

       % Acausal filtering of signal
       f1 = lower_cutoff; % low cutoff frequency
       f2 = upper_cutoff; % high cutoff frequency
        if f1 < 100 && f1 > 0
            tw = 10; % for ecog signal filtering
        else
            tw = 100; % transition width
        end
        
       att = 40; % stopband attenuation, in dB
       display = 'off'; % show bode plot?
       ht = bandpass(f1, f2, tw, fs, att, display); % FIR filter
       fftfiltfile(outfile{i}, 'temp1.raw', ht, 'int16');
       flipfile('temp1.raw', 'temp2.raw', 'int16');
       fftfiltfile('temp2.raw', 'temp3.raw', ht, 'int16');
       flipfile('temp3.raw', outfile{i}, 'int16');

       fclose('all');
       
   else
       fprintf('\n%s already exists\n', outfile{i});
   end

end % (for i)



trig_outfile = sprintf('%s-%s-%s-%s-%s-%s-fs%.0f-%s.raw', ...
   details.exp, details.site, details.depth(1:6), details.atten, details.stim,...
   details.probe, fs, trig_native_channel_name);% modified for simultaneous recording

d = dir(trig_outfile);

if ( isempty(d) )
    fprintf('%s\n%s\n\n', trigfile, trig_outfile);

%Reading and writing raw waveform data
fpathin = trigfile;
fpathout = trig_outfile; % one file change to use the addition

fidin = fopen(fpathin,'r');
fidout = fopen(fpathout,'w');

while ~feof(fidin)
   data = fread(fidin,nsamples,'uint16');
   if( ~isempty(data) )
      fwrite(fidout,data,'int16');
   end
end % (while)

fclose('all');

else
    fprintf('\n%s already exists\n', trig_outfile);
end
stim = details.stim;

return;






