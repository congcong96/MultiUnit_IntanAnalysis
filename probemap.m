%% for ASSY-156-ECoG-64B
output = [2:2:32 1:2:31 63:-2:33 64:-2:34]';
intan = [16:2:46 17:2:47 15:-2:1 63:-2:49 14:-2:0 62:-2:48]';
output2intan = [output, intan];
output2intan = sortrows(output2intan);

xposi = repmat([0 30 60 90],16,1);
yposi = repmat(115:30:(115+450),4,1);
yposi = yposi';

chnmap = [26 21 11 8;
35 34 64 61;
33 25 7 63;
24 38 60 10;
37 36 62 59;
28 23 9 6;
41 20 14 55;
40 46 52 58;
39 45 51 57;
27 32 2 5;
22 29 3 12;
44 18 16 54;
43 48 31 53;
42 50 1 56;
30 47 49 4;
19 17 15 13];

xposi = xposi(:);
yposi = yposi(:);
chnmap = chnmap(:);
chnmap = [xposi, yposi, chnmap];
chnmap = sortrows(chnmap,3);
chnmap = [chnmap output2intan(:,2)];