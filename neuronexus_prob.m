function [probinfo] = neuronexus_prob(probtype)
% neuronexus_prob(probtype) returns electrode positions of the Michigan 32 channel silicon probe.
%
% neuronexus_prob(probtype)
%
% probtype
%  Poly1 -> {'C5CA'} -> Model: A1x32-6mm-50-177-A32
%  Poly2 -> {'8B32', '8A5B'} -> Model: A1x32-Poly2-5mm-50s-177-A32
%  Poly3 -> {'638B', '638A'} -> Model: A1x32-Poly3-6mm-50s-177-A32
%
% updated to probtype from SerialNum and added Poly1 Natsumi 21 Jun 17
% Natsumi 27/02/17
% add ECoGx64B 
% Output:
%   probinfo.idxdepth: x position of the channel (lower left as origion)
%   probinfo.depth: distance from tip of electode
%   probinfo.posi_electrode: output channel number for the probe
%   probinfo.posi_intan: input channel number for intan headstage
% 2019-12-19, Congcong

%% check input and output
if  nargin ~= 1
    error('Input should be a type of Neuronexus probes (Choose from Poly1, Poly2, Poly3, LLNL or H31x64)');
end

% if (~ischar(nargin))
%     error('Input should be a string');
% end

%% pick up the model and get the configulation

% poly1 = {'C5CA'};
% poly2 = {'8B32', '8A5B'};
% poly3 = {'638B', '638A'};
%
% if ( sum( ismember(poly2, serialnum) ) )
%     model = 'A1x32-Poly2-5mm-50s-177-A32';
% elseif ( sum( ismember(poly3, serialnum) ) )
%     model = 'A1x32-Poly3-6mm-50s-177-A32';
% else
%     error('The serial number is not registerd');
% end

if ( sum( ismember([{'Poly1'},{'poly1'}], probtype) ) )
    model = 'A1x32-6mm-50-177-A32';
elseif ( sum( ismember([{'Poly2'},{'poly2'}], probtype) ) )
    model = 'A1x32-Poly2-5mm-50s-177-A32';
elseif ( sum( ismember([{'Poly3'},{'poly3'}, {'poly31x32'}], probtype) ) )
    model = 'A1x32-Poly3-6mm-50s-177-A32';
elseif ( sum( ismember([{'LLNL'},{'llnl'}], probtype) ) )
    model = 'Livermore-32Ch-Rat';
elseif ( sum( ismember([{'H3'},{'H31'},{'H31x64'}], probtype) ) )
    model = 'H31x64';
elseif (ismember({'ECoG64B'}, probtype))
    model = 'ASSY-156-ECoG-64B';
elseif sum(ismember( {'H22x32', 'H2'}, probtype))%H2 connected with neuronexus adaptor
    model = 'ASSY-77-H2';
elseif sum(ismember({'Tetrode1x64', 'Tetrode64'}, probtype))
    model = 'A4x4-tet';
elseif sum(ismember( {'smH22x32'}, probtype))%H2 connected with cambridge adaptor
    model = 'ASSY-77-H2-cambridge';
elseif sum(ismember({'poly21x48'}, probtype))
    model = 'A1x48-Poly2-5mm-100s-177';
elseif sum(ismember({'TetrodeB1x64'}, probtype))
    model = 'A4x4-tetB';
elseif sum(ismember({'Atlas1x16'}, probtype))
    model = 'Atlas1x16';
elseif sum(ismember({'smH2B2x32'}, probtype))
    model = 'ASSY-77-H2-cambridgeB';
else
    error('This type of prob is not registerd. (Choose from Poly1, Poly2, Poly3, LLNL or H31x64)');
end

[probmat] = getpositions(model);

%% extract the probe infor from probmat

% The distance from tip of electrode
posi_depth = probmat(:,1);
% The index for depth
posi_idxdepth = probmat(:,2);
% The amplifier/channel numbers for electrode
posi_electrode = probmat(:,3);
% The amplifier/channel numbers for intan
posi_intan = probmat(:,4);

if size(probmat) > 4
    posi_x = probmat(:,5);
else
    posi_x = zeros(size(probmat,1), 1);
end
%% save the probe info

probinfo.probtype = cell2mat(probtype);
probinfo.model = model;
probinfo.posi_idxdepth = posi_idxdepth;
probinfo.posi_depth = posi_depth;
probinfo.posi_intan = posi_intan;
probinfo.posi_electrode = posi_electrode;
probinfo.posi_x = posi_x;
end

%% [probmat] = getpositions(model);
function [probmat] = getpositions(model)

load('ElectrodePositions.mat')

if strcmp(model, 'A1x32-Poly2-5mm-50s-177-A32')
    IdxModel = 1;
elseif strcmp(model, 'A1x32-Poly3-6mm-50s-177-A32')
    IdxModel = 2;
elseif strcmp(model, 'A1x32-6mm-50-177-A32')
    IdxModel = 3;
elseif strcmp(model, 'Livermore-32Ch-Rat')
    IdxModel = 4;
elseif strcmp(model, 'H31x64')
    IdxModel = 5;
elseif strcmp(model, 'ASSY-156-ECoG-64B')
    IdxModel = 6;
elseif strcmp(model, 'ASSY-77-H2')
    IdxModel = 7;
elseif strcmp(model, 'A4x4-tet')
    IdxModel = 8;
elseif strcmp(model, 'ASSY-77-H2-cambridge')
    IdxModel = 9;
elseif strcmp(model, 'A1x48-Poly2-5mm-100s-177')
    IdxModel = 10;
elseif strcmp(model, 'A4x4-tetB')
    IdxModel = 11;
elseif strcmp(model, 'Atlas1x16')
    IdxModel = 12;
elseif strcmp(model,  'ASSY-77-H2-cambridgeB')
    IdxModel = 13;
end

probmat = ElectrodePositions(IdxModel).probmat;

end

%% to generate probmat
% 
% A32toOM32 = [32,30,31,28,29,27,25,22,23,21,26,24,20,19,18,17,1,4,13,14,15,16,11,9,7,5,3,2,6,8,10,12;32,31,30,29,28,27,26,25,24,23,17,18,19,20,21,22,16,15,14,13,12,11,1,2,3,4,5,6,7,8,9,10];
% OM32toIntan = [9,7,5,3,1,13,15,11,21,17,19,31,29,27,25,23,10,8,6,4,2,14,16,12,22,18,20,32,30,28,26,24;23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7];
% 
% % Poly1
% A32 = [17,16,18,15,19,14,20,13,21,12,22,11,23,10,24,9,25,8,26,7,27,6,28,5,29,4,30,3,31,2,32,1];
% % % Poly2
% % A32 = [10,9,8,7,6,5,4,3,2,1,11,12,13,14,15,16,23,24,25,26,27,28,29,30,31,32,22,21,20,19,18,17];
% % % Poly3
% % A32 = [3,2,1,4,5,6,7,8,9,10,11,17,16,18,15,19,14,20,13,21,12,30,31,32,29,28,27,26,25,24,23,22];
% 
% OM32 = NaN(size(A32)); Intan = NaN(size(A32));
% for idx = 1:length(A32)
%     
%     f = find(A32toOM32(1,:) == A32(idx));
%     OM32(idx) = A32toOM32(2,f);
%     f = find(OM32toIntan(1,:) == OM32(idx));
%     Intan(idx) = OM32toIntan(2,f);
%     
% end
% 
% Output = [A32;Intan]';
% ECoG = [26, 21, 11, 8,
%         35, 34, 64, 61,
%         33, 25, 7, 63,
%         24, 38, 60, 10,
%         37, 36, 62, 59,
%         28, 23, 9, 6,
%         41, 20, 14, 55,
%         40, 46, 52, 58,
%         39, 45, 51, 57,
%         27, 32, 2, 5,
%         22, 29, 3, 12,
%         44, 18, 16, 54,
%         43, 48, 31, 53,
%         42, 50, 1, 56,
%         30, 47, 49, 4,
%         19, 17, 15, 13]
% 
% ASSY-77-H2

% adaptorout = [34, 35, 62, 33, 60, 54, 57, 55, 10, 8, 11, 5, 32, 3, 30, 31,
%             64, 58, 63, 56, 61, 59, 52, 50, 15, 13, 6, 4, 9, 2, 7, 1,
%             53, 51, 49, 47, 45, 36, 37, 38, 27, 28, 29, 20, 18, 16, 14, 12,
%             48, 46, 44, 42, 40, 39, 43, 41, 24, 22, 26, 25, 23, 21, 19, 17];
% intanin = [16:2:46,
%             17:2:47,
%             15:-2:1, 63:-2:49,
%             14:-2:0, 62:-2:48];
% intanin = intanin';
% adaptorout = adaptorout';
% intan2adaptor = [intanin(:), adaptorout(:)];
% 
% probeout = [1:5, 7:2:15, 6:2:16, 59:-2:49, 64:-1:60, 58:-2:50;
%              17:21, 23:27, 22, 28, 32, 29, 30, 31, 43, 37, 33, 36, 35, 34, 48:-1:44, 42:-1:38];
% adaptorin = [24, 27, 22, 28, 29, 20, 18, 16, 14, 12, 26, 25, 23, 21, 19, 17, 39, 40, 42:2:48, 41, 38, 43, 37, 36, 45, 47:2:53;
%             10, 15, 8, 13, 6, 4, 9, 2, 7, 1, 11, 5, 32, 3, 30, 31, 54, 60, 33, 62, 35, 34, 55, 50, 57, 52, 59, 61, 56, 63, 58, 64];
% probeout = probeout';
% adaptorin = adaptorin';
% probe2adaptor = [adaptorin(:), probeout(:)];
% 
% intan2adaptor = sortrows(intan2adaptor, 2);
% probe2adaptor = sortrows(probe2adaptor, 1);
% probe2intan = [intan2adaptor, probe2adaptor(:,2)]; %[intan, adaptor, probe]
% 
% probe = [21, 23, 24, 30, 29, 16, 18, 20, 27, 19, 17, 25, 26, 32, 28, 22, 1, 3, 5, 7, 9, 11, 13, 15, 31, 14, 12, 10, 8, 6, 4, 2
%          64, 62, 60, 58, 56, 54, 52, 50, 34, 51, 53, 55, 57, 59, 61, 63, 44, 42, 41, 35, 36, 49, 47, 45, 38, 46, 48, 40, 39, 33, 37, 43];
% probe = probe';
% probe2position = [probe(:), 25*31-[0:25:25*31 0:25:25*31]', [zeros(32,1); 250*ones(32, 1)]];
% probe2intan = sortrows(probe2intan, 3);
% probe2position = sortrows(probe2position, 1);
% probe2intan = [probe2intan, probe2position(:,2:3)];
% probe2intan = sortrows(probe2intan, [5 4], {'ascend', 'descend'});
% probmat = [probe2intan(:,4),[1:64]' , probe2intan(:,3), probe2intan(:,1) probe2intan(:,5)];



%% A4x4-tet
% load('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
% ElectrodePositions(8).model = 'A4x4-tet';
% intanin = [30 31 28 29 26:-2:16   27:-2:17    37:2:47     32:36 38:2:46;
%            0:4 6:2:14      5:2:15    59:-2:49    62 63 60 61 58:-2:48];
% probeout = [41, 38, 43, 37, 36, 45, 47:2:53,       39, 40:2:48,    26, 25:-2:17,    24, 27, 22, 28, 29, 20, 18:-2:12;
%             55, 50, 57, 52, 59, 61, 56, 63, 58, 64,       54, 60, 33, 62, 35, 34,   11, 5, 32, 3, 30, 31,     10, 15, 8, 13, 6, 4, 9, 2, 7, 1];
% channel2intan = [probeout(:), intanin(:)];
% channel2intan = sortrows(channel2intan);
% probelayout = [9 12 8 6     10 13 7 4   11 15 5 2   14 16 3 1   
%     25 28 24 22     26 29 23 20     27 31 21 18     30 32 19 17     
%     41 44 40 38     42 45 39 36     43 47 37 34     46 48 35 33
%     57 60 56 54     58 61 55 52     59 63 53 50     62 64 51 49]';
% x = [[18 36 18 0     18 36 18 0      18 36 18 0  18 36 18 0]
%     [18 36 18 0     18 36 18 0      18 36 18 0  18 36 18 0]+200
%     [18 36 18 0     18 36 18 0      18 36 18 0  18 36 18 0]+400
%     [18 36 18 0     18 36 18 0      18 36 18 0  18 36 18 0]+600]';
% y = [[[536 518 500 518] [536 518 500 518]-150 [536 518 500 518]-300 [536 518 500 518]-450]
%     [[536 518 500 518] [536 518 500 518]-150 [536 518 500 518]-300 [536 518 500 518]-450]
%     [[536 518 500 518] [536 518 500 518]-150 [536 518 500 518]-300 [536 518 500 518]-450]
%     [[536 518 500 518] [536 518 500 518]-150 [536 518 500 518]-300 [536 518 500 518]-450]]';
% probeposition = [probelayout(:) x(:) y(:) [1:64]'];
% probeposition = sortrows(probeposition);
% probmap = [probeposition(:,3), probeposition(:,4), probeposition(:,1) channel2intan(:,2), probeposition(:,2)];
% probmap = sortrows(probmap, 2);
% ElectrodePositions(8).probmat = probmap;
% save('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')

 %% A1x48-Poly2
% load('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
% ElectrodePositions(10).model = 'A1x48-Poly2-5mm-100s-177';
% intanin = [30 31 28 29 26:-2:16   27:-2:17    37:2:47     32:36 38:2:46;
%            0:4 6:2:14      5:2:15    59:-2:49    62 63 60 61 58:-2:48];
% probeout = [25 0 27 0 0 29:2:37 0  0 26:2:32 0 0 23:-2:17 24 0 22 0 0 20:-2:12;
%     39 34 41 36 43 45 40 47 42 48 38 44 0 46 0 0 11 5 0 3 0 0 10 15 8 13 6 4 9 2 7 1];
% channel2intan = [probeout(:), intanin(:)];
% channel2intan = sortrows(channel2intan);
% probelayout = [10:-1:1 11:24 37:48 36:-1:25]';
% x = [zeros(24,1); 86*ones(24,1)];
% y = [[2300:-100:0]+75, [2300:-100:0]+75+50]';
% probeposition = [probelayout(:) x(:) y(:)];
% probeposition = sortrows(probeposition,3, 'descend');
% probeposition = [probeposition(:,1:3) [1:48]'];
% probeposition = sortrows(probeposition);
% probmap = [probeposition(:,3), probeposition(:,4), probeposition(:,1) channel2intan(17:end,2), probeposition(:,2)];
% probmap = sortrows(probmap, 2);
% ElectrodePositions(10).probmat = probmap;
% save('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
%% A4x4-tetB
%% ch16
% load('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
% ElectrodePositions(12).model = 'Atlas1x16';
% positiony = (0:150:2250)';
% positionx = zeros(16,1);
% probenum = (1:16)';
% intannum = [9 1 11 3 10 2 12 4 13 5 14 6 15 7 16 8]';
% probmap = [positiony, (1:16)', probenum, intannum, positionx];
% ElectrodePositions(12).probmat = probmap;
% save('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
%% smH22x32
% ASSY-77-H2
% load('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
% intanin = [ 32:36 38:2:46    37:2:47     27:-2:17          30 31 28 29 26:-2:16 ;
%            62 63 60 61 58:-2:48         59:-2:49      5:2:15      0:4 6:2:14  ];
% probeout = [1:5, 7:2:15, 6:2:16, 59:-2:49, 64:-1:60, 58:-2:50;
%              17:21, 23:27, 22, 28, 32, 29, 30, 31, 43, 37, 33, 36, 35, 34, 48:-1:44, 42:-1:38];
% probe2intan = [probeout(:), intanin(:)];
% probe2intan = sortrows(probe2intan);
% 
% probe = [21, 23, 24, 30, 29, 16, 18, 20, 27, 19, 17, 25, 26, 32, 28, 22, 1, 3, 5, 7, 9, 11, 13, 15, 31, 14, 12, 10, 8, 6, 4, 2
%          64, 62, 60, 58, 56, 54, 52, 50, 34, 51, 53, 55, 57, 59, 61, 63, 44, 42, 41, 35, 36, 49, 47, 45, 38, 46, 48, 40, 39, 33, 37, 43];
% probe = probe';
% probe2position = [probe(:), 25*31-[0:25:25*31 0:25:25*31]', [zeros(32,1); 250*ones(32, 1)], (1:64)'];
% probe2position = sortrows(probe2position);
% 
% probmat = [probe2position(:,2),probe2position(:,4) , probe2position(:,1), probe2intan(:,2), probe2position(:,3)];
% probmat = sortrows(probmat, 2);
% ElectrodePositions(9).probmat = probmat;
% save('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis/ElectrodePositions.mat', 'ElectrodePositions')
%% smH2B2x32
% ElectrodePositions(13).model =  'ASSY-77-H2-cambridgeB';
% probemat = ElectrodePositions(9).probmat;
% intan1 = [46:-2:16 47:-2:17 49:2:63 1:2:15 48:2:62 0:2:14];
% intan = [intan1; flip(intan1)]';
% intan = sortrows(intan);
% probemat = sortrows(probemat, 4);
% probemat(:,4) = intan(:,2);
% probemat = sortrows(probemat, 2);
% ElectrodePositions(13).probmat=  probemat;

%% to check wiring for LLNL probe
%
% A32 = [31 29 27 25 23 21 19 17 15 13 11 9 7 5 3 1 32 30 28 26 24 22 20 18 16 14 12 10 8 6 4 2];
% OM32s = [30 32 26 28 22 24 18 20 14 16 10 12 6 8 2 4; 29 31 25 27 21 23 17 19 13 15 9 11 5 7 1 3];
% Intans = [7 6 5 4 3 2 1 0 31 30 29 28 27 26 25 24; 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
%
% OM32toIntan = [OM32s(:) Intans(:)];
% 
% Intan = NaN(size(A32));
% for idx = 1:length(A32)
%     f = find(A32(idx) == OM32toIntan(:,1));
%     Intan(idx) = OM32toIntan(f,2);
% end
% 
% Output = [A32;Intan]';
