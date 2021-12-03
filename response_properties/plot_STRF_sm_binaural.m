nsite = 4;
site_idx = find([strf.site]==nsite);
mdb = 40;
dur = 20*60;
taxis = strf(1).taxis;
faxis = strf(1).faxis(152:256);
p = 0.002;
c = 1;
for ii = 1:length(site_idx)
    
    nplot = mod(ii,8);
    if nplot == 1
        f = figure('Renderer', 'painters', 'Position', [10 50 800 800]);
        sgtitle(sprintf('site%d', nsite))
    elseif nplot == 0
        nplot = 8;
    end
    rf_contra = strf(site_idx(ii)).rfcontra(152:256, :);
    rf_ipsi = strf(site_idx(ii)).rfipsi(152:256, :);
    subplot(4, 4, nplot*2-1)
    n0 = strf(site_idx(ii)).n0contra;
    plot_strf(rf_contra, n0, p, mdb, dur, taxis,faxis)
    title(sprintf('unit%d %dum',strf(site_idx(ii)).unit, strf(site_idx(ii)).position(2)))
    text(100, 50, sprintf('n0:%d', n0))
    subplot(4, 4, nplot*2)
    n0 = strf(site_idx(ii)).n0ipsi;
    plot_strf(rf_ipsi, n0, p, mdb, dur, taxis, faxis)
    
    if nplot == 8
        saveas(gcf, sprintf('binaural_sta_site%d_%d.jpg',nsite, c));
        c = c+1;
        close
    end
    
end
if nplot ~= 8
    saveas(gcf, sprintf('binaural_sta_site%d_%d.jpg',nsite, c));
    close
end