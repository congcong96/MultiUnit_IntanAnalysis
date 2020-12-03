function [BF, BW_max, BW_min, Latency,Q, tw, confirm, strf, faxis, taxis] = calc_bf_bw_latency(filter,taxis,faxis, varargin)
% calculate BF, BW, Latency for a filter (STA, MID1, MID")
% INPUTS
% filter: spectrotemporal receptive field
% taxis: time axis for the filter (ms)
% faxis: frequency axis for the filter (Hz)
% bw_crit: the threshold for obtaining the width of the distribution
% flag_ex: true or false, a flag for only using positive values. use with bw_crit = 0.9 (Atencio & Schreiner, 2008,  2013)
% OUTPUTS
% BF: best frequency (kHz)
% BW_max: maximum frequency of the bandwidth (kHz)
% BW_min: minimum frequency of the bandwidth (kHz)
% Latency: the peak in the time marginal (ms)
% BW: bandwidth (kHz), BW_max-BW_min
% Q: BF/BW
% Congcong, 2020-04-08
p = inputParser;
addParameter(p,'manual',0);
% to make the width of the distribution at 1/e of the peak height (or 0.75, Atencio et al. 2012)
addParameter(p,'bw_crit',0.37);
addParameter(p,'wsize',[1 1 200 150]);%100ms, 3 octave
addParameter(p,'rfsig',[]);%100ms, 3 octave

parse(p,varargin{:});
manual = p.Results.manual;
bw_crit = p.Results.bw_crit;
wsize = p.Results.wsize;
 
% peakstim = max(max(abs(filter(:, 101:301))));
% peakspon = max(max(abs(filter(:, 1:100))));
% %peakbefore spike is leass than 1.5 * peak after spike
% if peakspon > 0 && peakstim/peakspon < 1.5
%     [confirm, BF, BW_max, BW_min, Latency,Q, tw] = rejection(manual);
%     return
% end

faxislog = log2(faxis ./ min(faxis));  % change into log scale
faxislogu = linspace(min(faxislog),max(faxislog),1000); % up sample
if manual == 0
    [BF, fidx1, fidx2] = spectral_response(filter, [1, size(filter,2)], faxis, faxislog, faxislogu, bw_crit);
    
    %tVec = sum(fliplr(sta),1);
    fidx1d = round(fidx1/1000*size(filter,1));
    fidx2d = round(fidx2/1000*size(filter,1));
    if fidx1d == 316 || fidx2d == 1
        [confirm, BF, BW_max, BW_min, Latency,  Q, tw] = rejection(manual);
        return
    end
    if fidx1d == 0
        fidx1d = 1;
    end
    taxisu = linspace(min(taxis),max(taxis),1000);% upsample taxis
    [Latency, tidx1, tidx2]= temporal_response(filter, [fidx1d, fidx2d], taxis, taxisu, bw_crit);
    
    
    tidx1d = round(tidx1/1000*size(filter,2));
    tidx2d = round(tidx2/1000*size(filter,2));
    if tidx1d == 0
        tidx1d = 1;
    end
    [BF_tmp, fidx1_tmp, fidx2_tmp] = spectral_response(filter, [tidx1d tidx2d], faxis, faxislog, faxislogu, bw_crit);
    c = 1;
    while BF_tmp ~= BF && c < 10
        BF = BF_tmp;
        fidx1 = fidx1_tmp;
        fidx2 = fidx2_tmp;
        fidx1d = round(fidx1/1000*size(filter,1));
        fidx2d = round(fidx2/1000*size(filter,1));
        if fidx1d == 0
            fidx1d = 1;
        end
        [Latency, tidx1, tidx2]= temporal_response(filter, [fidx1d fidx2d], taxis, taxisu, bw_crit);
        tidx1d = round(tidx1/1000*size(filter,2));
        tidx2d = round(tidx2/1000*size(filter,2));
        if tidx1d == 0
            tidx1d = 1;
        end
        [BF_tmp, fidx1_tmp, fidx2_tmp] = spectral_response(filter, [tidx1d tidx2d], faxis, faxislog, faxislogu, bw_crit);
        c = c + 1;
    end
    BW_max = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
    BW_min = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
    if BW_max/BW_min < 2^0.095 %noise when bw is smaller than half modulation cycle
        [confirm, BF, BW_max, BW_min, Latency,  Q, tw] = rejection(manual);
        return
    end
    tw =  [taxisu(tidx1), taxisu(tidx2)];
    if Latency <= 0 || Latency >= 100 || diff(tw) < 0.0038
        [confirm, BF, BW_max, BW_min, Latency,  Q, tw] = rejection(manual);
        return
    end
    
    Q = BF/ (BW_max-BW_min);
    
elseif manual
    rfsig = p.Results.rfsig;
    subplot(122)
    plot_strf_raw(rfsig, faxis, taxis)
    subplot(121)
    plot_strf_raw(filter, faxis, taxis)
    [BF, fidx1, fidx2] = spectral_response(filter, [1, size(filter,2)], faxis, faxislog, faxislogu, bw_crit);
    fidx1d = round(fidx1/1000*size(filter,1));
    fidx2d = round(fidx2/1000*size(filter,1));
    if fidx1d == 316 || fidx2d == 1
        [confirm, BF, BW_max, BW_min, Latency,  Q, tw] = rejection(manual);
        if confirm
            return
        end
    end
    if fidx1d == 0
        fidx1d = 1;
    end
    taxisu = linspace(min(taxis),max(taxis),1000);% upsample taxis
    [Latency, tidx1, tidx2, idxLatency]= temporal_response(filter, [fidx1d, fidx2d], taxis, taxisu, bw_crit);
    if Latency <= 0 || Latency > 100
        [confirm, BF, BW_max, BW_min, Latency,  Q, tw] = rejection(manual);
        if confirm
            return
        end
    end
    tidx1d = round(tidx1/1000*size(filter,2));
    tidx2d = round(tidx2/1000*size(filter,2));
    idxLatency = round(idxLatency/1000*size(filter,2));
    if tidx1d == 0
        tidx1d = 1;
    end
    [BF_tmp, fidx1_tmp, fidx2_tmp, idxBF] = spectral_response(filter, [tidx1d tidx2d], faxis, faxislog, faxislogu, bw_crit);    
    fidx1d = round(fidx1_tmp/1000*size(filter,1));
    fidx2d = round(fidx2_tmp/1000*size(filter,1));
    idxBF = round(idxBF/1000*size(filter,1));
    if rfsig(idxLatency, idxBF) == 0
        [confirm, BF, BW_max, BW_min, Latency,  Q, tw] = rejection(manual);
        if confirm
            return
        end
    end
    hold on
    plot([1 idxLatency], [idxBF idxBF], 'k-');%draw the ling for BF
    plot([idxLatency idxLatency], [1 idxBF],'k-');
    plot([1 size(filter,2)], [fidx1d fidx1d], 'k--');%lower bound of BW
    plot([1 size(filter,2)], [fidx2d fidx2d], 'k--');%upper bounf for BW
    plot([tidx1d tidx1d], [1 size(filter,1)], 'k--');%lower bound of BW
    plot([tidx2d tidx2d], [1 size(filter,1)], 'k--');%upper bounf for BW
    
    roi = images.roi.Rectangle(gca,'Position', wsize);
    idxLatency = round(idxLatency/1000*size(filter,2));
    idxBF = round(idxBF/1000*size(filter,1));
    
   
   
    confirm = input('\nAcceptance confirm? 1=yes, 0=no: ');
else
    confirm = NaN;
end
end

function [BF, idx1, idx2, idxBF]= spectral_response(filter, flim, faxis,faxislog, faxislogu, bw_crit)
    fVec = sum(abs(filter(:, flim(1):flim(2))), 2);
    fVecu = interp1(faxislog,fVec,faxislogu,'spline');
    [M, idxBF] = max(fVecu);
    BF = (min(faxis)*2^faxislogu(idxBF))/1000; %kHz
    p1 = 1;
    p2 = 1;
    while idxBF > p1 && abs(fVecu(idxBF - p1)) > bw_crit*M
        p1 = p1+1;
    end
    while idxBF + p2 < length(fVecu) && abs(fVecu(idxBF + p2)) > bw_crit*M
        p2 = p2+1;
    end
    idx1 = idxBF - p1 + 1;
    idx2 = idxBF + p2 - 1;
end

function [Latency, idx1, idx2, idxLatency]= temporal_response(filter, tlim, taxis, taxisu, bw_crit)
    tVec = sum(filter(tlim(1):tlim(2),:));
    tVecu = interp1(taxis,tVec,taxisu,'spline');
    [M, idxLatency] = max(abs(tVecu));
    Latency = taxisu(idxLatency)*1000; %ms
    p1 = 1;
    p2 = 1;
    while idxLatency > p1 && abs(tVecu(idxLatency - p1)) > bw_crit*M
        p1 = p1+1;
    end
    while idxLatency + p2 < length(tVecu) && abs(tVecu(idxLatency + p2)) > bw_crit*M
        p2 = p2+1;
    end
    idx1 = idxLatency - p1 + 1;
    idx2 = idxLatency + p2 - 1;
end

function  [confirm, BF, BW_max, BW_min, Latency, Q, tw] = rejection(manual)
    if manual 
        confirm = input('\nrejection confirm? 1=yes, 0=no: ');
    else
        confirm = NaN;
    end
    BF = NaN;
    BW_max = NaN;
    BW_min = NaN;
    Latency = NaN;
    Q = NaN;
    tw = NaN;
end