sessions = dir('E:\Congcong\Documents\emsemble_thalamus\*CH');

savefolder = 'E:\Congcong\Documents\emsemble_thalamus\figure\multiunit\summary';
stimfolder = 'E:\Congcong\Documents\stimulus\thalamus';
fra10stim = 'freq_resp_area_stimulus_flo500Hz_fhi32000Hz_nfreq21_natten8_nreps10_fs96000_param';
fra10params = load(fullfile(stimfolder,fra10stim));

%% FRA10 properties (used later for MGB subnucleus categorization)
<<<<<<< HEAD
for ii = 1:length(sessions)
=======
for ii = 18:length(sessions)
>>>>>>> 3919e44a1bc045d3a781b69b9806350b2623d932
    sessionfolder = fullfile(sessions(ii).folder, sessions(ii).name);
    cd(sessionfolder)
    
    fra10recordings = dir('site*fra10*');
    for jj = 1:length(fra10recordings)
        % load matrix
        fra10folder = fullfile(fra10recordings(jj).folder,fra10recordings(jj).name);
        fra10file = dir(fullfile(fra10folder, '*thresh-fra10.mat'));
        load(fullfile(fra10file.folder, fra10file.name), 'thresh', 'trigger')
        
        % FRA properties
        cd(fra10folder)
        rasterfile = dir('*-raster.mat');
        load(rasterfile.name, 'raster')
        
        % order raster to match the order of channels in thresh
        raster = OrderChannel(thresh, raster);
        
        
        fraproperties = [];
        for kk = 1:length(thresh)
            [tctmp, frapropertiestmp, psth, delay, rpw, n0, HW] = fra_properties(thresh(kk), raster(kk), trigger, fra10params);
            thresh(kk).psth = psth;
            thresh(kk).delay = delay;
            thresh(kk).tc = tctmp;
            thresh(kk).rpw = rpw;
            thresh(kk).HW = HW;
            thresh(kk).n0 = n0;
            
            frapropertiestmp.probe = thresh(kk).probe;
            frapropertiestmp.chan = thresh(kk).chan;
            fraproperties = [fraproperties frapropertiestmp];
        end
        
        save(fra10file.name, 'thresh', 'fraproperties', 'raster', '-append')
        % plot fra and psth
        %plot_fra_psth(thresh, fraproperties, fullfile(savefolder,[fra10file.name(1:27) 'fra_tc_psth']), 1)
    end
end
%% dmrrep properties
% RI is the correlation of 2 halves of PSTH, 100 repeat
% 1000 repeat for null distribution
for ii = 1:length(sessions)
    sessionfolder = fullfile(sessions(ii).folder, sessions(ii).name);
    cd(sessionfolder)
    
    dmrreprecordings = dir('site*_dmrrep_*');
    for jj = 1:length(dmrreprecordings)
        % load matrix
        dmrrepfolder = fullfile(dmrreprecordings(jj).folder,dmrreprecordings(jj).name);
        dmrrepfile = dir(fullfile(dmrrepfolder, '*thresh-raster.mat'));
        load(fullfile(dmrrepfile.folder,dmrrepfile.name), 'thresh', 'raster')
        fprintf('Processing %s\n', dmrrepfile.name)
        
        % order thresh and raster
        raster = OrderChannel(thresh, raster);
        
        % check Reliability of raster
        dmrrepRI = batch_calc_dmrrep_RI(raster);
        save(fullfile(dmrrepfile.folder, dmrrepfile.name), 'raster', 'dmrrepRI', '-append');
    end
end
%% dmr properties
% reliability index calculated on banshee
p = 0.002;
mdb = 40;
<<<<<<< HEAD
for ii = 1:length(sessions)
=======
for ii = 1%:length(sessions)
>>>>>>> 3919e44a1bc045d3a781b69b9806350b2623d932
    sessionfolder = fullfile(sessions(ii).folder, sessions(ii).name);
    cd(sessionfolder)
    
    dmrrecordings = dir('site*_dmr_*');
    for jj = 1:length(dmrrecordings)
        % load matrix
        dmrfolder = fullfile(dmrrecordings(jj).folder,dmrrecordings(jj).name);
        dmrfile = dir(fullfile(dmrfolder, '*thresh-strf.mat'));
        load(fullfile(dmrfile.folder,dmrfile.name), 'thresh', 'trigger', 'strf')
        dur = (trigger(end)-trigger(1))/20000;
        % order thresh and strf
        if isfield(strf, 'probe')
            probe = {thresh.probe};
            strf.probe = probe{:};
        end
        strf = OrderChannel(thresh, strf);
        
        % get properties of STRF features
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
        strfproperties = [];
        
        
        for kk = 1:length(strf)
            rf = strf(kk).rfcontra;
            n0 = strf(kk).n0contra;
            rfsig = significant_strf(rf, p, n0, mdb, dur);
<<<<<<< HEAD
            rfsig = rfsig;
            [BF, BW_max, BW_min, Latency, Q, tw, confirm] = calc_bf_bw_latency(rfsig,taxis,faxis);
=======
            [BF, BW_max, BW_min, Latency, Q, tw] = calc_bf_bw_latency(rfsig,taxis,faxis, 'manual', 1, 'rfraw', rf);
>>>>>>> 3919e44a1bc045d3a781b69b9806350b2623d932
            strfproperties(kk).chan = strf(kk).chan;
            %strfproperties(kk).probe = strf(kk).probe;
            strfproperties(kk).BF = BF;
            strfproperties(kk).Latency = Latency;
            strfproperties(kk).BW = [BW_min, BW_max];
            strfproperties(kk).Q = Q;
            strfproperties(kk).tw = tw;
        end
        
        if isfield(strf, 'probe')
            probe = {strf.probe};
            strfproperties.probe = probe{:};
        end
        
        rtf = strf_parameters(strf, dur);

        save(fullfile(dmrfile.folder, dmrfile.name),'strf', 'rtf', 'strfproperties', '-append');

    end
end

%% plor MU properties
p = 0.002;
mdb = 40;
for ii = 1:length(sessions)
    sessionfolder = fullfile(sessions(ii).folder, sessions(ii).name);
    cd(sessionfolder)
    
    fra10recordings = dir('site*fra10*');
    dmrrecordings = dir('site*_dmr_*');
    dmrreprecordings = dir('site*dmrrep*');
    for jj = 1:length(fra10recordings)
        
        fra10folder = fullfile(fra10recordings(jj).folder,fra10recordings(jj).name);
        fra10file = dir(fullfile(fra10folder, '*thresh-fra10.mat'));
        load(fullfile(fra10file.folder, fra10file.name), 'thresh', 'fraproperties')
        
        dmrfolder = fullfile(dmrrecordings(jj).folder,dmrrecordings(jj).name);
        dmrfile = dir(fullfile(dmrfolder, '*thresh-strf.mat'));
        load(fullfile(dmrfile.folder,dmrfile.name), 'strf', 'rtf', 'strfproperties', 'trigger')
        
        dmrrepfolder = fullfile(dmrreprecordings(jj).folder,dmrreprecordings(jj).name);
        dmrrepfile = dir(fullfile(dmrrepfolder, '*thresh-raster.mat'));
        load(fullfile(dmrrepfile.folder,dmrrepfile.name), 'dmrrepRI', 'raster')
        
        probe = {thresh.probe};
        [strf.probe] = probe{:};
        [rtf.probe] = probe{:};
        [strfproperties.probe] = probe{:};
        dur = (trigger(end)-trigger(1))/20000;
        strf = OrderChannel(thresh, strf);
        rtf = OrderChannel(thresh, rtf);
        strfproperties = OrderChannel(thresh, strfproperties);
        dmrrepRI = OrderChannel(thresh,dmrrepRI);
        raster = OrderChannel(thresh, raster);
        
        fignum = 1;
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
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
            if isnan(fraproperties(kk).Q10)
                title(sprintf('%dum      %dspk', thresh(kk).position(2), thresh(kk).n0))
            else
                hold on
                plot(fraproperties(kk).sboundry, 'k', 'linewidth', 2)
                title(sprintf('%dum  %.1fkHz  %dspk', thresh(kk).position(2), fraproperties(kk).CF, thresh(kk).n0))
            end
            
            %plot PSTH
            subplot(4, 6, 6*nplot - 4)
            bar(thresh(kk).psth, 'k', 'EdgeColor', 'k')
            hold on
            rpw = thresh(kk).rpw;
            bar(rpw(1)+1:rpw(2), thresh(kk).psth(rpw(1)+1:rpw(2)), 'b', 'EdgeColor', 'b')
            if ~isnan(thresh(kk).delay)
                title(sprintf('%d:%s%d %dms', kk, thresh(kk).amplifier, thresh(kk).chan, thresh(kk).delay))
            else
                title(sprintf('%d:%s%d', kk, thresh(kk).amplifier, thresh(kk).chan))
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
            
            if any(~isnan(strfproperties(kk).BF))
                hold on
                for mm = 1:length(strfproperties(kk).BF)
                    if isnan(strfproperties(kk).BF(mm))
                        continue
                    end
                    latency = strfproperties(kk).Latency(mm);
                    idxLatency = find(taxis>latency/1000, 1);
                    BF =  strfproperties(kk).BF(mm);
                    idxBF = find(faxis>=BF*1000, 1);
                    BW = strfproperties(kk).BW([mm*2-1, mm*2]);
                    fidx2d =  find(faxis>=BW(2)*1000, 1);
                    fidx1d = find(faxis>=BW(1)*1000, 1);
                    
                    plot([1 idxLatency], [idxBF idxBF], 'k-');%draw the line for BF
                    plot([idxLatency idxLatency], [1 idxBF],'k-');
                    if isempty(fidx2d) || isempty(fidx1d)
                        continue
                    end
                    plot([1 size(rfsig,2)], [fidx1d fidx1d], 'k--');%lower bound of BW
                    plot([1 size(rfsig,2)], [fidx2d fidx2d], 'k--');%upper boun for BW
                    
                   
                    
                    %title(sprintf('%.0fms %.1fkHz [%.1f, %.1f]',latency, BF, BW(1), BW(2)))
                end
            end
            
            
            %plot rtf
            
            if ~isnan(strfproperties(kk).BF)
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
                rtf(kk).rtfparams = RTFparam;
                plot(RTFparam.tBMF, RTFparam.sBMF,'k+','linewidth',2)
                axis([min(Fm) max(Fm) 0 max(RD)])
                title(sprintf('sBMF %.1f/tBMF %.f',RTFparam.sBMF, RTFparam.tBMF ))
                %     axis([0 max(Fm)/2 0 max(RD)/2])
            end
            %plot response to dmrrep
            subplot(4, 6, [6*nplot-1, 6*nplot])
            plotSpikeRaster(raster(kk).raster_mat>0,'PlotType','vertline');
            set(gca,'xtick', 0:500:5000, 'xticklabel',0:10);
            title(sprintf('RI %.2f    p=%.3f', dmrrepRI(kk).RI, dmrrepRI(kk).p))
            
            if nplot == 4
                saveas(gca, fullfile(savefolder, sprintf('%s-site%d-%s-%d.jpg', thresh(1).exp, thresh(1).site, thresh(1).probe, fignum)))
                close
                fignum = fignum + 1;
            end
        end
        save(fullfile(dmrfile.folder, dmrfile.name),'rtf', '-append')
    end
end
