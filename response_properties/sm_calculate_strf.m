function strf = calculate_strf_v2(spk, trigger, envfile, tbefore, tafter)
% procedures same as calculate STRF, but preload stim_mat as global
% variable
global stim_mat
if nargin == 3
    T1 = 0.050; % 50 ms
    T2 = 0.200; % 200 ms
else
    T1 = tafter;
    T2 = tbefore;
end

i2 = strfind(envfile,'SM');
i3 = strfind(envfile,'TM');
i4 = strfind(envfile,'db');

sm = str2double(envfile(i2(end)-1));
tm = str2double(envfile(i3(end)-2:i3(end)-1));
mdb = str2double(envfile(i4(end)-2:i4(end)-1));

try
    atten = spk.atten;
catch
    atten = 0;
end

fs = double(spk(1).fs);
spl = 105 - atten;

%Loading Parameter Data
index=strfind(envfile,'.spr');
ParamFile=[envfile(1:index(1)-1) '_param.mat'];
load(ParamFile, 'Fs', 'DF', 'NF', 'faxis')

RMSP=-mdb/2;					% RMS value of normalized Spectral Profile
PP=mdb^2/8;					% Modulation Depth Variance

for i = 1:length(spk)
    fprintf('Processing spk%d/%d ......\n', i, length(spk))
    sp = spk(i).spiketimes;
    spet = round(sp / 1000 * double(fs)); % convert ms to sample number
    
    %----------------------calculate STRF----------------------------
    
    %Converting Temporal Delays to Sample Numbers
    N1=round(T1*Fs/DF); % DF = temporal downsampling factor for spectral profile
    N2=round(T2*Fs/DF);
    
    %Align spikes to onsett of stimulus
    onset  = trigger(1);
    spet = spet - onset;
    spet = floor(spet* Fs / fs /DF );
    spet(spet < N2) = [];
    spet(spet > length(stim_mat)-N1) = [];
    
    
    
    STRF = zeros(NF, N1+N2);
    for k=1:length(spet)
        STRF=STRF+ stim_mat(:, spet(k)-N2+1:spet(k)+N1)*mdb - RMSP;
    end
    
    n0 = length(spet);
    T = (length(stim_mat) - N1 - N2) / Fs * DF ;
    w0 = n0/T;
    taxis=(-N1:N2-1)/(Fs/DF);
    STRF=w0/PP*fliplr(STRF)/n0;
    %-------------------fill in strf structure-------------------------
    strf(i).chan = spk(i).chan;
    
    if isfield(spk, 'position')
        strf(i).position = spk(i).position;
    end
    strf(i).sm = sm;
    strf(i).tm = tm;
    strf(i).mdb = mdb;
    strf(i).spl = spl;
    strf(i).taxis = taxis;
    strf(i).faxis = faxis;
    strf(i).n0contra = n0;
    strf(i).w0contra = w0;
    strf(i).rfcontra = STRF; % use single precision to save memor
    strf(i).fs = fs;
end




