function [strct2] = OrderChannel(strct1, strct2)
% strct2 is ordered based on order of channels in strct1

if length(strct1) == 128 && any([strct1.chan] ~= [strct2.chan])
    f = cellfun(@(x) strcmp(x, strct1(1).probe), {strct2.probe});
    strct2tmp = strct2(f);
    chanidx = [1:64; [strct1(1:64).chan]]';
    chanidx = sortrows(chanidx, 2);
    strct2tmpidx =[1:64; [strct2tmp(1:64).chan]]';
    strct2tmpidx = sortrows(strct2tmpidx, 2);
    strct2tmpidx = [chanidx(:,1), strct2tmpidx(:,1)];
    strct2tmpidx = sortrows(strct2tmpidx, 1);
    strct2tmp = strct2tmp(strct2tmpidx(:,2));
    
    f = cellfun(@(x) strcmp(x, strct1(65).probe), {strct2.probe});
    strct2(65:end) = strct2(f);
    chanidx = [65:128; [strct1(65:end).chan]]';
    chanidx = sortrows(chanidx, 2);
    strct2idx =[65:128; [strct2(65:end).chan]]';
    strct2idx = sortrows(strct2idx, 2);
    strct2idx = [chanidx(:,1), strct2idx(:,1)];
    strct2idx = sortrows(strct2idx, 1);
    strct2(65:end) = strct2(strct2idx(:,2));
    
    strct2(1:64) = strct2tmp;

elseif length(strct1) == 64 && any([strct1.chan] ~= [strct2.chan])
    chanidx = [1:64; strct1.chan]';
    chanidx = sortrows(chanidx, 2);
    strct2idx = [1:64; strct2.chan]';
    strct2idx = sortrows(strct2idx, 2);
    strct2idx = [chanidx(:,1), strct2idx(:,1)];
    strct2idx = sortrows(strct2idx, 1);
    strct2 = strct2(strct2idx(:,2));
end