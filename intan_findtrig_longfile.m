function [trigger, triptrig] = intan_findtrig_longfile(trigpath, Thresh)
%findtrig_longfile  Find triggers from a .raw or .mich int16 file
%
%  [trigger, triptrig] = findtrig_longfile(trigpath, Thresh, M)
%  This function is used to find the triggers within a long trigger file,
%  where multiple stimuli were presented. This function differs from 
%  findtrig.m in that it returns the position of the triple triggers, and
%  does not remove the extra triggers for each triple trigger.
%
%  trigpath : filename 
%
%  Thresh : Threshhold : [-1 , 1]. Default: Thresh==.75
%
%  M : Block Size. Default: M==1024*256. Unless you really care, I don't
%     recommend supplying this input argument. The default works nicely.
%
%  trigger : Returned Trigger Time Vector (in sample number)
%  triptrig : vector of triple trigger times (in sample number)
%
% Craig Atencio
% 6-8-2012
%
% [trigger, triptrig] = findtrig_longfile(trigpath, Thresh, M);


%Checking Inputs
if nargin == 1
	Thresh = .75;
end 


% fileinfo = dir('board-ADC-00.dat');
% num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
% fid = fopen('board-ADC-00.dat', 'r');
% v = fread(fid, num_samples, 'uint16');
% fclose(fid);
% v = v * 0.000050354; % convert to volts


%Checking if trigpath is a string
if ( ~ischar(trigpath) )
    error('TRIGPATH must be a string.');
else
    %Remove leading and trailing whitespaces
    trigpath = strtrim(trigpath);
end

fileinfo = dir(trigpath);
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(trigpath, 'r');
X = fread(fid, num_samples, 'uint16');
fclose(fid);
MeanX = mean(X);
MaxX = max(X);
X = (X-MeanX)/MaxX;
ThreshN = max(X)*Thresh; %max(X)*0.70; % add new threshold 

% Finding Triggers
trigger = [];

X( X >= ThreshN ) = 1;
X( X < ThreshN ) = 0;
D = diff(X);
trigger = find( D >= abs(Thresh) ) + 1;	%Finding Trigger Locations
trigger = trigger';
%Re-calculating inter-trigger spacing
Dtrigger = diff(trigger);
maxD = median(Dtrigger); %get trigger interval

%Finding triple trigger(s)
triple_ind = find(Dtrigger(1:end-1) < 0.34*maxD & Dtrigger(2:end) < 0.34*maxD );
%triple_ind = find(Dtrigger(1:end-2) < 0.34*maxD & Dtrigger(2:end-1) < 0.34*maxD & Dtrigger(3:end) < 0.34*maxD); % chenged natsumi 10Oct17
triptrig = trigger(triple_ind); % this should tell us where the triple triggers are


%Removing additional triggers if extra triple triggers are found
if(isempty(triple_ind))

    disp('Could not find initial triple trigger.');

else %Triple triggers found

   % Multiple triple triggers found. We will assume that the first and 
   % last triple trigger were for the same stimulus. This is the case
   % for the static ripple/moving ripple stimuli, but not for other stimuli
   if diff(triple_ind) > 1 %( length(triple_ind) == 2 )  % changed  this to avoid over trimming when triple_ind = [1,2] Natsumi
     trigger = trigger( triple_ind(1):triple_ind(2) );

   %There was only one triple trigger
   else
     trigger = trigger(triple_ind(1):end);
   end
   
   if ( length(triple_ind) == 2 &&diff(triple_ind) > 1) % add this instead of above lines Natsumi
       trigger = [trigger(1) trigger(5:end)];
   else
       %Making triple trigger count as one trigger
       trigger = [trigger(1) trigger(4:end)];
   end

end

return;



