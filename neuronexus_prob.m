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
elseif ( sum( ismember([{'Poly3'},{'poly3'}], probtype) ) )
    model = 'A1x32-Poly3-6mm-50s-177-A32';
elseif ( sum( ismember([{'LLNL'},{'llnl'}], probtype) ) )
    model = 'Livermore-32Ch-Rat';
elseif ( sum( ismember([{'H3'},{'H31'},{'H31x64'}], probtype) ) )
    model = 'H31x64';
elseif (ismember({'ECoG64B'}, probtype))
    model = 'ASSY-156-ECoG-64B';
elseif sum(ismember([{'H2'}, {'H22x32'}], probtype))
    model = 'ASSY-77-H2';
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
