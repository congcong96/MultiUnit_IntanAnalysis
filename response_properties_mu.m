sessions = dir('E:\Congcong\Documents\emsemble_thalamus\*CH');
savefolder = 'E:\Congcong\Documents\emsemble_thalamus\figure\multiunit\summary';

stimfolder = 'E:\Congcong\Documents\stimulus\thalamus';
fra10stim = 'freq_resp_area_stimulus_flo500Hz_fhi32000Hz_nfreq21_natten8_nreps10_fs96000_param';
fra10params = load(fullfile(stimfolder,fra10stim));

for ii = 1:length(sessions)
    sessionfolder = fullfile(sessions(ii).folder, sessions(ii).name);
    cd(sessionfolder)
    
    fra10recordings = dir('site*fra10*');
    for jj = 1:length(fra10recordings)
        %% load matrix
        cd(sessionfolder)
        fra10folder = fullfile(fra10recordings(jj).folder,fra10recordings(jj).name);
        fra10file = dir(fullfile(fra10folder, '*thresh-fra10.mat'));
        load(fullfile(fra10file.folder, fra10file.name), 'thresh', 'trigger')
        site = thresh(1).site;
        depth = thresh(1).depth;
        
        dmrfolder = dir(sprintf('site%d_%dum*_dmr_*', site, depth));
        dmrfile = dir(fullfile(dmrfolder.folder, dmrfolder.name, '*-strf.mat'));
        load(fullfile(dmrfile.folder, dmrfile.name), 'strf')
        
        dmrrepfolder = dir(sprintf('site%d_%dum*_dmrrep_*', site, depth));
        dmrrepfile = dir(fullfile(dmrrepfolder.folder, dmrrepfolder.name, '*-raster.mat'));
        load(fullfile(dmrrepfile.folder, dmrrepfile.name), 'raster')
        dmrraster = raster;
        
        %% FRA properties
        cd(fra10folder)
        rasterfile = dir('*-raster.mat');
        load(rasterfile.name, 'raster')
        
        % order raster to mach the order of channels in thresh
        if length(raster) == 128 && raster(1).chan == 0
            f = cellfun(@(x) strcmp(x, thresh(1).probe), {raster.probe});
            rastertmp = raster(f);
            strftmp = strf(f);
            dmrrastertmp = dmrraster(f);
            chanidx = [thresh(1:64).chan]+1;
            rastertmp = rastertmp(chanidx);
            strftmp = strftmp(chanidx);
            probe = {rastertmp.probe};
            strftmp.probe = probe{:};
            dmrrastertmp = dmrrastertmp(chanidx);
            f = cellfun(@(x) strcmp(x, thresh(65).probe), {raster.probe});
            raster(65:end) = raster(f);
            probe = {raster(65:end).probe};
            strf(65:end) = strf(f);
            strf(65:end).probe = probe{:};
            dmrraster(65:end) = dmrraster(f);
            chanidx = [thresh(65:end).chan]+65;
            raster(65:end) = raster(chanidx);
            strf(65:end) = strf(chanidx);
            dmrraster(65:end) = dmrraster(chanidx);
            raster(1:64) = rastertmp;
            strf(1:64) = strftmp;
            dmrraster(1:64) = dmrrastertmp;
        elseif raster(1).chan == 0
            chanidx = [thresh.chan]+1;
            raster = raster(chanidx);
            strf = strf(chanidx);
            dmrraster = dmrraster(chanidx);
        end
        
        fraproperties = [];
        for kk = 1:length(thresh)
            [tctmp, frapropertiestmp, psth, delay, rpw, n0] = fra_properties(thresh(kk), raster(kk), trigger, fra10params);
            thresh(kk).psth = psth;
            thresh(kk).delay = delay;
            thresh(kk).tc = tctmp;
            thresh(kk).rpw = rpw;
            thresh(kk).n0 = n0;
            
            frapropertiestmp.probe = thresh(kk).probe;
            frapropertiestmp.chan = thresh(kk).chan;
            fraproperties = [fraproperties frapropertiestmp];
        end
        
        save(fra10file.name, 'thresh', 'fraproperties', 'raster', '-append')
        
        %% strf properties
        dur = 60*15; 
        rtf = strf_parameters(strf, dur);
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
        chan = {strf.chan};
        strfproperties = [];
        for kk = 1:length(strf)
            [BF, BW_max, BW_min, Latency, Q, confirm] = calc_bf_bw_latency(strf(kk).rfcontra,taxis,faxis);
            strfproperties(kk).BF = BF;
            strfproperties(kk).Latency = Latency;
            strfproperties(kk).BW = [BW_min, BW_max];
            strfproperties(kk).Q = Q;
            strfproperties(kk).confirm = confirm;
            
        end
        [strfproperties.chan] = chan{:};
        if isfield(strf, 'probe')
            probe = {strf.probe};
            strfproperties.probe = probe{:};
        end
        save(fullfile(dmrfile.folder, dmrfile.name),'strf', 'rtf', 'strfproperties', '-append');
        %% dmrrep
        dmrrepRI = batch_calc_dmrrep_RI(dmrraster);
        raster = dmrraster;
        save(fullfile(dmrrepfile.folder, dmrrepfile.name), 'raster', 'dmrrepRI', '-append');
        
        %% plor MU properties
        nfig = 0;
        p = 0.002;
        mdb = 40;
        
        rtfparams = [];
        fignum = 1;
        for kk = 1:length(fraproperties)
            nplot = mod(kk ,4);
            if nplot == 1
                figure('Renderer', 'painters', 'Position', [30 30 1500 1000]);
            elseif nplot == 0
                nplot = 4;
            end
            
            %plot FRA
            subplot(4, 6, 6*nplot - 5);
            tuning_curve_plot(thresh(kk).tc);
            if isnan(fraproperties(kk).CF)
                title(sprintf('%dum      %dspk', thresh(kk).position(2), thresh(kk).n0))
            else
                plot(fraproperties(kk).boundry, 'w', 'linewidth', 3)
                title(sprintf('%dum  %dkHz  %dspk', thresh(kk).position(2), fraproperties(kk).CF, thresh(kk).n0))
            end
            
            %plot PSTH
             subplot(4, 6, 6*nplot - 4)
             bar(thresh(kk).psth, 'k', 'EdgeColor', 'k')
             hold on
             rpw = thresh(kk).rpw;
             bar(thresh(kk).psth(rpw(1)+1:rpw(2)), 'b', 'EdgeColor', 'b')
             if ~isnan(thresh(kk).delay)
                title(sprintf('%dms', thresh(kk).delay))
             end
             
             %plot strf
             subplot(4, 6, 6*nplot - 3)
             rf = strf(kk).rfcontra;       
             n0 = strf(kk).n0contra;
             w0 = strf(kk).w0contra;
             [rfsig] = significant_strf(rf, p, n0, mdb, dur);
             h = fspecial('gaussian',3,3);
             rfsig = imfilter(rfsig, h);
             imagesc(rfsig)
             cmap = cschemes('rdbu', 21);
             colormap(gca, cmap)
             clim = max(abs(rfsig(:)));
             caxis([-clim, clim]);
             axis xy
             xticks([1:100:500 500])
             xticklabels(-50:50:200)
             xlabel('ms')
             yticks(1:50:316)
             yticklabels([0.5 1 2 4 8 16 32])
             ylabel('kHz')
             
             if ~isnan(strfproperties(kk).BF)
                 hold on
                 latency = strfproperties(kk).Latency;
                 idxLatency = find(taxis>latency/1000, 1);
                 BF =  strfproperties(kk).BF;
                 idxBF = find(faxis>BF*1000, 1);
                 BW = strfproperties(kk).BW;
                 fidx1d = find(faxis>BW(1)*1000, 1);
                 fidx2d =  find(faxis>BW(2)*1000, 1);
                 
                 plot([1 idxLatency], [idxBF idxBF], 'k-');%draw the ling for BF
                 plot([1 size(rfsig,2)], [fidx1d fidx1d], 'k--');%lower bound of BW
                 plot([1 size(rfsig,2)], [fidx2d fidx2d], 'k--');%upper bounf for BW
                 plot([idxLatency idxLatency], [1 idxBF],'k-');
                 
                 title(sprintf('%.0fms %.1fkHz [%.1f, %.1f]',latency, BF, BW(1), BW(2)))
             end
            
             
             %plot rtf    
             subplot(4, 6, 6*nplot - 2)
             idx_energy = 4;
             cmap = cschemes('rdbu', 21);
             RTF = rtf(kk).rtf(:,:,idx_energy);
             Fm = rtf(kk).tmf;
             RD = rtf(kk).xmf;
             
             Max=max(max(RTF));
             imagesc(Fm,RD,RTF,[0 Max]),shading flat,colormap(gca, cmap) %axis square
             set(gca,'Ydir','normal')
             hold on
             contour(Fm,RD,RTF,[.25 .75]*Max,'k');
             
             [RTFparam]=rtf_parameters04(RTF,RD,Fm);
             rtfparams = [rtfparams, RTFparam];
             plot(RTFparam.tBMF, RTFparam.sBMF,'k+','linewidth',2)
             axis([min(Fm) max(Fm) 0 max(RD)])
             title(sprintf('sBMF %.1f/tBMF %.f',RTFparam.sBMF, RTFparam.tBMF ))
             %     axis([0 max(Fm)/2 0 max(RD)/2])
             
             %plot response to dmrrep
             subplot(4, 6, [6*nplot-1, 6*nplot])
             plotSpikeRaster(raster(kk).raster_mat>0,'PlotType','vertline');
             set(gca,'xtick', 0:500:5000, 'xticklabel',0:10);
             title(sprintf('RI %.2f    p=%.3f', dmrrepRI(kk).RI, dmrrepRI(kk).p))
             
             if nplot == 4
                 saveas(gca, fullfile(savefolder, sprintf('%s-site%d-%s-%d.jpg', thresh(1).exp, thresh(1).site, thresh(1).probe, fignum)))
                 fignum = fignum+1;
                 close 
             end
        end
    end
    save(fullfile(dmrfile.folder, dmrfile.name),'rtfparams', '-append')
end