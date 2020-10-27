function dmrrepRI = batch_calc_dmrrep_RI(raster)

dmrrepRI = struct('exp',{raster.exp}, 'chan',{raster.chan}, 'probe', {raster.probe});

rep = 100;
permmat = zeros(rep, 50);
for ii = 1:rep
    permmat(ii,:) = randperm(50);
end

for ii = 1:length(raster)
    rastermat = raster(ii).raster_mat;
    [RI, p, corr_raster, corr_null] = calc_dmrrep_SI(rastermat, permmat);
    dmrrepRI(ii).RI = RI;
    dmrrepRI(ii).p = p;
    dmrrepRI(ii).corr_raster = corr_raster;
    dmrrepRI(ii).corr_null = corr_null;
end
end

function  [RI, p, corr_raster, corr_null] = calc_dmrrep_SI(rastermat, permmat)
    
    permmatA = permmat;
    permmatA(permmatA<=25) = 1;
    permmatA(permmatA>25) = 0;
    rasterA = permmatA*rastermat;
    
    permmatB = permmat;
    permmatB(permmatB<=25) = 0;
    permmatB(permmatB>25) = 1;
    rasterB = permmatB*rastermat;
    
    corr_raster = corr(rasterA', rasterB');
    corr_raster = diag(corr_raster);
    RI = mean(corr_raster);
    
    rep = 1000;
    corr_null = zeros(rep, 100);
    for ii = 1:rep
        if mod(ii, 50) == 0
            fprintf('iteratio %d/1000... \n', ii)
        end
        shift = floor(rand(50,1)*length(rastermat));
        rastermatnull =rastermat;
        for jj = 1:size(rastermat,1)
            rastermatnull(jj,:) = circshift(rastermat(jj,:), shift(jj));
        end
        rasternullA = permmatA*rastermatnull;
        rasternullB = permmatB*rastermatnull;
        corr_null(ii,:) = diag(corr(rasternullA', rasternullB'));
    end
    RI_null = mean(corr_null, 2);
    p = sum(RI_null >= RI)/rep;

end