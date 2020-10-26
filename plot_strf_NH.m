function plot_strf_NH(strf, trigger, probinfo, depth, savepath, flim, tlim)
%plot_strf - Plots STRF data obtained using the Michigan 16
%   channel silicon probe.
%
%   plot_strf(strf, trigger, flim, tlim)
%
%   strf is a struct array holding the receptive field data
%
%   trigger is a vector of trigger times associated with the
%     ripple stimlui
%
%   flim and tlim are optional input arguments.
%
%   flim is an optional 1x2 vector, holding the frequency range
%      over which the strfs are plotted. If flim is not input then
%      the strfs are plotted over the full frequency range.
%
%   tlim is the same as flim, except it is over the time axis.
%
%   caa 6/28/02
%
% modified plot_strf(strf, trigger, flim, tlim) little bit, just to change some vidualization.
% Natsumi 240217
% deleted exp, depth, site ...
% Congcong 11/07/2019
% modiefied to print out pdf as specified by savepath (folder/name)
% Congcong 11222019
if ( nargin < 0 || nargin > 6 )
    error('You need 2 to 4 input args.');
end

if ( nargin == 5 )
    flim = [];
    tlim = [];
end

if ( nargin == 6 )
    tlim = [];
end

time = 1000 .* strf(1).taxis; % change to ms
freq = strf(1).faxis;


% Now get some time tick marks
if ( isempty(tlim) )
    t = round(min(time)):abs(round(min(time))):max(time);
    ttick = [];
    for i = round(min(time)):abs(round(min(time))):max(time)
        temp = abs(i+0.01-time);
        ttick = [ttick find(temp == min(temp))];
    end
else
    
    dttick = ( abs(tlim(1)) + abs(tlim(2)) ) / 5;
    t = round(tlim(1)):dttick:round(tlim(2));
    
    index_tlim_min = find( time >= tlim(1) );
    index_tlim_min = min(index_tlim_min);
    
    index_tlim_max = find( time <= tlim(2) );
    index_tlim_max = max(index_tlim_max);
    
    ttick = [];
    for i = round(tlim(1)):dttick:round(tlim(2))
        temp = abs(i+0.01-time);
        ttick = [ttick find(temp == min(temp))];
    end
    ttick = ttick - (index_tlim_min - 1);
end


% Now we get some frequency tick marks
if isempty(flim)
    nocts = floor(log2(max(freq)/min(freq)));
    f = round(min(freq).*2.^(0:2:nocts));
    ftick = [];
    for i = 0:2:nocts
        temp = abs(i+0.001-log2(freq/min(freq)));
        ftick = [ftick find(temp == min(temp))];
    end
else
    
    nocts = floor(log2(max(flim(2))/min(flim(1))));
    f = round(min(flim(1)).*2.^(0:1:nocts));
    ftick = [];
    for i = 0:1:nocts
        temp = abs(i+0.001-log2(freq/min(flim(1))));
        ftick = [ftick find(temp == min(temp))];
    end
    index_flim_min = find( freq >= flim(1) );
    index_flim_max = find( freq <= flim(2) );
    index_flim_min = min(index_flim_min);
    index_flim_max = max(index_flim_max);
    ftick = ftick - (index_flim_min - 1);
end


fs = strf(1).fs; % sampling rate of A/D system
mdb = strf(1).mdb;
sigma2 = (mdb)^2 / 8; % variance of dynamic moving ripple
sigma = sqrt(sigma2); % std of dynamic moving ripple
dur = ( trigger(end)-trigger(1) ) / fs; % total duration of ripple, in sec

probidx = probinfo.posi_idxdepth;
probidxintan = probinfo.posi_intan;
len = length(strf);
fignum = 1;

gaussian = fspecial('gaussian', [20 20], 7);
% imagesc(gaussian)
% pause
figurec = 1;
for i = 1:len
    
    if mod(i,16) == 1 %i <16
        fignum = 1;
        figure('Renderer', 'painters', 'Position', [30 30 1000 1000]);
    else
        fignum = fignum + 1;
    end
    f2 = find(i == probidx);
    chan = strf(probidxintan(f2)+1).chan;
    
    
    deep = depth-probinfo.posi_depth(f2);
    
    subplot(4,4,fignum);
    
    rf = double(strf(probidxintan(f2)+1).rfcontra);
    
    p = 0.002;
    n0 = strf(probidxintan(f2)+1).n0contra;
    w0 = strf(probidxintan(f2)+1).w0contra;
    [rfsig] = significant_strf(rf, p, n0, mdb, dur);
    
    if ( ~isempty(flim) )
        rfsig = rfsig(index_flim_min:index_flim_max, :);
    end
    
    if ( ~isempty(tlim) )
        rfsig = rfsig(:, index_tlim_min:index_tlim_max);
    end
    
    %    rfsig = imfilter(rfsig,gaussian);
    
    if isempty(rfsig)
        rfsig = zeros(length(freq), length(time));
    end
  
    minmin = min(min(rfsig));
    maxmax = max(max(rfsig));
    boundary = max([abs(minmin) abs(maxmax)]);
    if  maxmax == 0
        BF=NaN;
    else
        
        [~, midx] = max(sum(rfsig,1)); %find time bin when the total energy is largest;
        
        
        BF = freq(rfsig(:,midx) == max(rfsig(:,midx)));
    end

    
    imagesc(rfsig);
%     cmap = cschemes('rdbu', 21);
%     colormap(cmap);
    colormap('jet');
    set(gca,'ydir', 'normal');
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    
    if  fignum == 1
        set(gca,'xtick', ttick, 'xticklabel',t);
        xlabel('Time (ms)');
    else
        set(gca,'xtick', ttick, 'xticklabel','');
    end
    
    if  fignum == 1
        set(gca,'ytick', ftick, 'yticklabel',f/1000);
        ylabel('Freq (kHz)');
    else
        set(gca,'ytick', ftick, 'yticklabel','');
    end
    
    %title(sprintf('%.0f:%.0f-%.0f %.0fum BF %.1f kHzn_0%.0f w_0%.2f', i, chan, model, deep, BF/1000, n0, w0));
    title(sprintf('%.0f:A%.0f %.0fum', i, chan, deep),'FontName','Arial','FontSize',9);
    if mod(i,16) == 0
        saveas(gcf, [savepath sprintf('-%d.png', figurec)]);
        figurec = figurec + 1;
        close
    end
end

end
