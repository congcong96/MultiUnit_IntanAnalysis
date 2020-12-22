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