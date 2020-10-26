%TUNING_CURVE_PLOT plot tuning curve data
%   TUNING_CURVE_PLOT(TC) plots the tuning curve data stored in the tuning
%   curve struct TC (see tuning_curve_calculate for fields).
%
%   Written by Jonathan Shih 2008-03-11
%   add BF (largest response collapsed on frequency axis)
%   only values cross threshold is colored
%   Congcong
function bf = tuning_curve_plot(tc, threshold)

%Checking input data
if(length(tc) > 1)
    error('Can only plot one tuning curve at a time');
end

if nargin < 2
    threshold = 0;
end

freq = tc.freq;
%Plotting tuning curve
newmat = smoothmat(tc.tcmat);
idx = newmat < threshold;
newmat(idx) = 0;
imagesc(newmat);
axis xy
xticks(round(linspace(1,length(freq),4)));
xticklabels(round(freq(round(linspace(1,length(freq),4)))/1000,1));
if length(tc.atten) == 15
    yticks(5:5:15);
    yticklabels(75+tc.atten(5:5:15));
else
    yticks(2:3:8);
    yticklabels(75+tc.atten(2:3:8));
end

tc_f = sum(newmat);
[~, I] = max(tc_f);
bf = tc.freq(I);
%cmap = cbrewer('seq', 'Blues', 9);
%colormap(gca, cmap);
colormap(gca, 'jet')

if length(freq) == 21
    text(12, 2,sprintf('%.1fk',bf/1000),'FontSize',10,'FontWeight','bold', 'Color', 'w')
else
    text(20, 2,sprintf('%.1fk',bf/1000),'FontSize',10,'FontWeight','bold', 'Color', 'w')
end
%tickpref;
% xlabel('Frequency (kHz)');
% ylabel('Attenuation (dB)');
