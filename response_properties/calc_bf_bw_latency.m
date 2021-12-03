function [BF, BW_max, BW_min, Latency,Q, tw, sig] = calc_bf_bw_latency(filter,taxis,faxis, varargin)
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

% mannual mode instruction:
% when a rf is rejected because the frequency or latency is out of range or
% the maximum feature is too small to be significant, it need to be
% confirmed by keyboard input
% when the rf is not rejected, keyboard input is also required to proceed
% then user select ROI for further analysis
%
% Congcong, 2020-04-08
p = inputParser;
% to make the width of the distribution at 1/e of the peak height (or 0.75, Atencio et al. 2012)
addParameter(p,'bw_crit',0.37);
addParameter(p,'manual',0);
addParameter(p,'rfraw',[]);
parse(p,varargin{:});
bw_crit = p.Results.bw_crit;
manual =  p.Results.manual;
rfraw =  p.Results.rfraw;

faxislog = log2(faxis ./ min(faxis));  % change into log scale
faxislogu = linspace(min(faxislog),max(faxislog),1000); % up sample
taxisu = linspace(min(taxis),max(taxis),1000);% upsample taxis
if max(abs(filter(:))) < 3*std(filter(:))
    rejection
    sig = 0;
    return
end
if manual
    subplot(131)
    plot_strf_raw(rfraw, faxis, taxis)
    subplot(132)
    plot_strf_raw(filter, faxis, taxis)
    [reject, BF, BW_max, BW_min, Latency,  Q, tw, sig] = auto_rejection(filter,faxis,taxis, bw_crit);
    if reject
        confirm =  input('\nInsignificant STRF, proceed? 1=yes(analyze anyway), 0 = no (discard): ');
        if isempty(confirm)|| confirm == 0
            return
        end
    end
    
    pfilter = filter;
    nfilter = filter;
    pfilter(filter<0) = 0;
    nfilter(filter>0) = 0;
    %% for positive response
    
    subplot(133)
    hold off
    plot_strf_raw(pfilter, faxis, taxis)
    [reject, BFe, BW_maxe, BW_mine, Latencye,  Qe, twe, sige] = auto_rejection(pfilter,faxis,taxis, bw_crit);
    if  reject
        confirm =  input('\nInsignificant STRF, proceed? 1=yes(analyze anyway), 0=no(discard): ');
        if confirm
            [BFe, fidx1, fidx2, Latencye, tidx1, tidx2] = find_peak(pfilter, faxis, taxis,  bw_crit);
            BW_maxe = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
            BW_mine = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
            
            twe =  [taxisu(tidx1), taxisu(tidx2)];
            Qe = BFe/ (BW_maxe-BW_mine);
        end
    else
        [BFe, fidx1, fidx2, Latencye, tidx1, tidx2] = find_peak(pfilter, faxis, taxis,  bw_crit);
        BW_maxe = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
        BW_mine = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
        
        twe =  [taxisu(tidx1), taxisu(tidx2)];
        Qe = BFe/ (BW_maxe-BW_mine);
    end
    if ~isnan(BFe)
        sige = 1;
        draw_BF_Latency(pfilter, taxis, faxis, BFe,BW_mine, BW_maxe, Latencye, twe)
        
        confirm =  input('\nNext step? 0=discard, 1=accept, 2=rescale:');
        
        while confirm == 2
            subplot(132)
            roi = drawrectangle(gca);
            position = round(roi.Position);
            if position(1)<=0
                position(1) = 1;
            end
            ptaxis = taxis(position(1): position(1)+position(3)-1);
            pfaxis = faxis(position(2): position(2)+position(4)-1);
            pfaxislog = log2(pfaxis ./ min(pfaxis));  % change into log scale
            pfaxislogu = linspace(min(pfaxislog),max(pfaxislog),1000); % up sample
            ptaxisu = linspace(min(ptaxis),max(ptaxis),1000);% upsample taxis
            filter = pfilter(position(2): position(2)+position(4)-1, position(1): position(1)+position(3)-1);
            
            [BFe, fidx1, fidx2, Latencye, tidx1, tidx2] = find_peak(filter, pfaxis, ptaxis,  bw_crit);
            BW_maxe = (min(pfaxis)*2^pfaxislogu(fidx2))/1000; %kHz
            BW_mine = (min(pfaxis)*2^pfaxislogu(fidx1))/1000; %kHz
            
            twe =  [ptaxisu(tidx1), ptaxisu(tidx2)];
            Qe = BFe/ (BW_maxe-BW_mine);
            subplot(133)
            plot_strf_raw(filter, pfaxis, ptaxis)
            draw_BF_Latency(filter, ptaxis, pfaxis, BFe,BW_mine, BW_maxe, Latencye, twe)
            delete(roi)
            confirm =  input('\nNext step? 0=discard, 1=accept, 2=rescale:');
        end
        
        if isempty(confirm) || confirm == 0
            [~, BFe, BW_maxe, BW_mine, Latencye,  Qe, twe] = rejection;
        end
    end
    %% for negtive response
    subplot(133)
    hold off
    plot_strf_raw(nfilter, faxis, taxis)
    [reject,  BFn, BW_maxn, BW_minn, Latencyn,  Qn, twn, sign] = auto_rejection(nfilter,faxis,taxis, bw_crit);
    if  reject
        confirm =  input('\nInsignificant STRF, proceed? 1=yes(analyze anyway), 0=no(discard): ');
        if confirm
            [BFn, fidx1, fidx2, Latencyn, tidx1, tidx2] = find_peak(nfilter, faxis, taxis,  bw_crit);
            BW_maxn = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
            BW_minn = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
            
            twn =  [taxisu(tidx1), taxisu(tidx2)];
            Qn = BFn/ (BW_maxn-BW_minn);
        end
    else
        [BFn, fidx1, fidx2, Latencyn, tidx1, tidx2] = find_peak(nfilter, faxis, taxis,  bw_crit);
        BW_maxn = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
        BW_minn = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
        
        twn =  [taxisu(tidx1), taxisu(tidx2)];
        Qn = BFn/ (BW_maxn-BW_minn);
    end
    
    if  ~isnan(BFn)
        sign = 1;
        draw_BF_Latency(nfilter, taxis, faxis, BFn,BW_minn, BW_maxn, Latencyn, twn)
        
        
        confirm =  input('\nNext step? 0=discard, 1=accept, 2=rescale:');
        
        while confirm == 2
            
            subplot(132)
            roi = drawrectangle(gca);
            position = round(roi.Position);
            if position(1)<=0
                position(1) = 1;
            end
            ntaxis = taxis(position(1): position(1)+position(3)-1);
            nfaxis = faxis(position(2): position(2)+position(4)-1);
            nfaxislog = log2(nfaxis ./ min(nfaxis));  % change into log scale
            nfaxislogu = linspace(min(nfaxislog),max(nfaxislog),1000); % up sample
            ntaxisu = linspace(min(ntaxis),max(ntaxis),1000);% upsample taxis
            filter = nfilter(position(2): position(2)+position(4)-1, position(1): position(1)+position(3)-1);
            
            [BFn, fidx1, fidx2, Latencyn, tidx1, tidx2] = find_peak(filter, nfaxis, ntaxis,  bw_crit);
            BW_maxn = (min(nfaxis)*2^nfaxislogu(fidx2))/1000; %kHz
            BW_minn = (min(nfaxis)*2^nfaxislogu(fidx1))/1000; %kHz
            
            twn =  [ntaxisu(tidx1), ntaxisu(tidx2)];
            Qn = BFn/ (BW_maxn-BW_minn);
            subplot(133)
            plot_strf_raw(filter, nfaxis, ntaxis)
            draw_BF_Latency(filter, ntaxis, nfaxis, BFn,BW_minn, BW_maxn, Latencyn, twn)
            delete(roi)
            confirm =  input('\nNext step? 0=discard, 1=accept, 2=rescale:');
        end
        
        if isempty(confirm)|| confirm == 0
            [~, BFn, BW_maxn, BW_minn, Latencyn,  Qn, twn,sign] = rejection;
        end
    end
    
    BF = [BFe, BFn];
    BW_max = [BW_maxe; BW_maxn];
    BW_min = [BW_mine; BW_minn];
    Latency = [Latencye, Latencyn];
    Q = [Qe, Qn];
    tw = [twe; twn];
    sig = [sige, sign];
else
    [reject, BF, BW_max, BW_min, Latency,  Q, tw, sig] = auto_rejection(filter,faxis,taxis, bw_crit);
    if reject
        return
    end
    
    pfilter = filter;
    nfilter = filter;
    pfilter(filter<0) = 0;
    nfilter(filter>0) = 0;
    
    
    [reject, BFe, BW_maxe, BW_mine, Latencye,  Qe, twe, sige] = auto_rejection(pfilter,faxis,taxis, bw_crit);
    if  ~reject
        [BFe, fidx1, fidx2, Latencye, tidx1, tidx2] = find_peak(pfilter, faxis, taxis,  bw_crit);
        BW_maxe = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
        BW_mine = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
        
        twe =  [taxisu(tidx1), taxisu(tidx2)];
        Qe = BFe/ (BW_maxe-BW_mine);
    end
    
    [reject, BFn, BW_maxn, BW_minn, Latencyn,  Qn, twn, sign] = auto_rejection(nfilter,faxis,taxis, bw_crit);
    if  ~reject
        [BFn, fidx1, fidx2, Latencyn, tidx1, tidx2] = find_peak(nfilter, faxis, taxis,  bw_crit);
        BW_maxn = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
        BW_minn = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
        twn =  [taxisu(tidx1), taxisu(tidx2)];
        Qn = BFn/ (BW_maxn-BW_minn);
    end
    
    
    BF = [BFe, BFn];
    BW_max = [BW_maxe; BW_maxn];
    BW_min = [BW_mine; BW_minn];
    Latency = [Latencye, Latencyn];
    Q = [Qe, Qn];
    tw = [twe; twn];
    sig = [sige, sign];
end

end

function [BF, idx1, idx2, idxBF]= spectral_response(filter, flim, faxis,faxislog, faxislogu, bw_crit)
    fVec = sum(abs(filter(:, flim(1):flim(2))), 2);
    fVecu = interp1(faxislog,fVec,faxislogu,'spline');
    [M, idxBF] = max(fVecu);
    BF = (min(faxis)*2^faxislogu(idxBF))/1000; %kHz
    p1 = 1;
    p2 = 1;
    while idxBF > p1+1 && abs(fVecu(idxBF - p1)) > bw_crit*M && abs(fVecu(idxBF - p1)) >  abs(fVecu(idxBF - p1-1)) 
        p1 = p1+1;
    end
    while idxBF + p2+1 < length(fVecu) && abs(fVecu(idxBF + p2)) > bw_crit*M  &&  abs(fVecu(idxBF + p2)) > abs(fVecu(idxBF + p2+1))
        p2 = p2+1;
    end
    idx1 = idxBF - p1 + 1;
    idx2 = idxBF + p2 - 1;
end

function [Latency, idx1, idx2, idxLatency]= temporal_response(filter, tlim, taxis, taxisu, bw_crit)
    tVec = sum(filter(tlim(1):tlim(2),:),1);
    tVecu = interp1(taxis,tVec,taxisu,'spline');
    [M, idxLatency] = max(abs(tVecu));
    Latency = taxisu(idxLatency)*1000; %ms
    p1 = 1;
    p2 = 1;
    while idxLatency > p1+1 && abs(tVecu(idxLatency - p1)) > bw_crit*M %&&  abs(tVecu(idxLatency - p1)) >  abs(tVecu(idxLatency - p1-1))
        p1 = p1+1;
    end
    while idxLatency + p2+1 < length(tVecu) && abs(tVecu(idxLatency + p2)) > bw_crit*M %&& abs(tVecu(idxLatency + p2)) > abs(tVecu(idxLatency + p2+1))
        p2 = p2+1;
    end
    idx1 = idxLatency - p1 + 1;
    idx2 = idxLatency + p2 - 1;
end

function  [reject, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection
    reject = 1;
    BF = NaN;
    BW_max = NaN;
    BW_min = NaN;
    Latency = NaN;
    Q = NaN;
    tw = [NaN NaN];
    sig = 0;
end

function  [reject, BF, BW_max, BW_min, Latency,  Q, tw, sig] = auto_rejection(filter,faxis,taxis, bw_crit)
faxislog = log2(faxis ./ min(faxis));  % change into log scale
faxislogu = linspace(min(faxislog),max(faxislog),1000); % up sample
[~, fidx1, fidx2] = spectral_response(filter, [1, size(filter,2)], faxis, faxislog, faxislogu, bw_crit);
fidx1d = round(fidx1/1000*size(filter,1));
fidx2d = round(fidx2/1000*size(filter,1));
if fidx2d >= length(faxis)-1 || fidx1d <= 1%if the best frequency is out of bound, reject
    %fprintf('\n Frequency out of bound\n')
    [reject, BF, BW_max, BW_min, Latency,  Q, tw, sig] = rejection;
    if reject
        return
    end
end

taxisu = linspace(min(taxis),max(taxis),1000);% upsample taxis
[Latency, tidx1, tidx2]= temporal_response(filter, [fidx1d, fidx2d], taxis, taxisu, bw_crit);
if Latency <= 0 || Latency > 100%if latency is smaller than 0 or larger than 100ms, reject
    fprintf('latency out of range\n')
    [reject, BF, BW_max, BW_min, Latency,  Q, tw, sig] = rejection;
    if reject
        return
    end
end

tidx1d = round(tidx1/1000*size(filter,2));
tidx2d = round(tidx2/1000*size(filter,2));
if tidx1d == 0
    tidx1d = 1;
end
[BF, fidx1, fidx2] = spectral_response(filter, [tidx1d tidx2d], faxis, faxislog, faxislogu, bw_crit);
tw =  [taxisu(tidx1), taxisu(tidx2)];
BW_max = (min(faxis)*2^faxislogu(fidx2))/1000; %kHz
BW_min = (min(faxis)*2^faxislogu(fidx1))/1000; %kHz
if diff(tw) < 0.0026 || BW_max/BW_min < 2^0.095 %noise when bw is smaller than half modulation cycle
    fprintf('feature too small\n')
    [reject, BF, BW_max, BW_min, Latency,  Q, tw, sig] = rejection;
    return
end
Q = BF/ (BW_max-BW_min);
reject = 0;
sig = 1;
end

function [BF, fidx1, fidx2, Latency, tidx1, tidx2] = find_peak(filter, faxis, taxis,  bw_crit)
faxislog = log2(faxis ./ min(faxis));  % change into log scale
faxislogu = linspace(min(faxislog),max(faxislog),1000); % up sample
[BF, fidx1, fidx2] = spectral_response(filter, [1, size(filter,2)], faxis, faxislog, faxislogu, bw_crit);
fidx1d = round(fidx1/1000*size(filter,1));
fidx2d = round(fidx2/1000*size(filter,1));

taxisu = linspace(min(taxis),max(taxis),1000);% upsample taxis
[Latency, tidx1, tidx2]= temporal_response(filter, [fidx1d, fidx2d], taxis, taxisu, bw_crit);
tidx1d = round(tidx1/1000*size(filter,2));
if tidx1d == 0
    tidx1d = 1;
end
tidx2d = round(tidx2/1000*size(filter,2));
    
[BF_tmp, fidx1_tmp, fidx2_tmp] = spectral_response(filter, [tidx1d tidx2d], faxis, faxislog, faxislogu, bw_crit);
c = 0;
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
end

function draw_BF_Latency(filter, taxis, faxis, BF,BW_min, BW_max, Latency, tw)
hold on
idxLatency = find(taxis>Latency/1000,1);
tidx1d = find(taxis>=tw(1),1);
tidx2d = find(taxis>=tw(2),1);
idxBF = find(faxis>=BF*1000,1);
fidx1d =  find(faxis>=BW_min*1000,1);
fidx2d =  find(faxis>=BW_max*1000,1);
if isempty( fidx2d )
    fidx2d = size(filter,1);
end
if isempty(idxLatency)
    idxLatency = size(filter,2);
end
plot([1 idxLatency], [idxBF idxBF], 'k-');%draw the line for BF
plot([idxLatency idxLatency], [1 idxBF],'k-');
plot([1 size(filter,2)], [fidx1d fidx1d], 'k--');%lower bound of BW
plot([1 size(filter,2)], [fidx2d fidx2d], 'k--');%upper bounf for BW
plot([tidx1d tidx1d], [1 size(filter,1)], 'k--');%lower bound of BW
plot([tidx2d tidx2d], [1 size(filter,1)], 'k--');%upper bounf for BW
hold off
end
