function raster = plot_raster_mu(raster, probinfo, depth, savepath)

probidx = probinfo.posi_idxdepth;
probidxintan = probinfo.posi_intan;

fignum = 1;
xt = 0:10;
for i = 1:length(raster)
    if mod(i,4) == 1 %i <16
        fignum = 1;
        figure('Renderer', 'painters', 'Position', [30 30 1000 1000]);
    else
        fignum = fignum + 1;
    end
    
    f2 = i == probidx;
    chan = raster(probidxintan(f2)+1).chan;
    deep = depth-probinfo.posi_depth(f2);
    raster(chan+1).depth = deep;
    
    subplot(4,1,fignum);
    plotSpikeRaster(raster(chan+1).raster_mat>0,'PlotType','vertline');
    title(sprintf('%.0f:A%.0f %.0fum', i, chan, deep),'FontName','Arial','FontSize',9);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    if mod(i,4) == 0
        set(gca,'xtick', 0:500:5000, 'xticklabel',xt);
        xlabel('time (s)');
        ylabel('trial')
        export_fig(savepath,'-pdf', '-append');
        close
    end
end
end