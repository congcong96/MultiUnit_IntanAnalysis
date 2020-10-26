function [trigger, outfile] = findtrig(trigpath, Thresh, M)
%findtrig  Find triggers from a .raw or .mich int16 file
%
%  [trigger, outfile] = findtrig(trigpath, Thresh, M)
%
%  trigpath : filename 
%
%  Thresh : Threshhold : [-1 , 1]. Default: Thresh==.75
%
%  M : Block Size. Default: M==1024*256. Unless you really care, I don't
%     recommend supplying this input argument. The default works nicely.
%
%  trigger : Returned Trigger Time Vector (in sample number)
%
%  outfile : output file name. The trigger file name without the .raw or
%            .mich ending. This makes it easy to save the trigger data
%            without having to copy and paste a file name.
%
% Craig Atencio
% 6-8-2012

narginchk(1,3)

%Checking Inputs
if nargin<2
	M = 1024*256;
	Thresh = .6;
elseif ( nargin < 3 )
	M = 1024*256;
end 



%Checking if trigpath is a string
if ( ~isstr(trigpath) )
    error('TRIGPATH must be a string.');
else
    %Remove leading and trailing whitespaces
    trigpath = strtrim(trigpath);
end

%If trigpath is a .ncs file, converts it to a .raw file
if ( lower(trigpath((end-2):end)) == 'ncs' )
    ncs2raw(trigpath,[trigpath(1:(end-3)) 'raw'],0);
    trigpath = [trigpath(1:(end-3)) 'raw'];
end

outfile = regexprep(trigpath, '\.raw', '\.mat');

% Opening trigger file
fid = fopen(trigpath,'r');

% Finding the Aproximate Mean in the File. Opening file and finding the 
% Mean by Searching 10 Segments
MeanX=0;
for j=1:10
   %Reading Data
   X = fread(fid,M,'int16')';
   MeanX = MeanX + mean(X);
end
MeanX = MeanX / 10;


%Finding the Max in the File. 
frewind(fid);

MaxX = 0;
while ~feof(fid)
   %Reading Data
   X = fread(fid,M,'int16')';
   MaxX = max([MaxX abs(X-MeanX)]);
end
fclose(fid);


%Closing Files
fclose('all');


% Finding Triggers
trigger = [];
count = 0;

%Opening trigger file
fid = fopen(trigpath,'r');

%Loading data and finding Triggers
while ~feof(fid)

   % Reading Input File - normalize to between [-1, 1]
	X = ( fread(fid, M,'int16')'-MeanX ) / MaxX;

   if ( Thresh > 0 )
      X(X<0) = 0;
   else
      X(X>0) = 0;
      X = -X;
   end


   % Set the triggers to a maximum/minimum value to eliminate the Gibbs
   % phenomenon effect
   X( X > 0.5 ) = 1;

% clf;
% plot(X)
% pause


   %Adding First Element to avoid missing trigger 
   %if trigger falls exactly on edge
   if ( count == 0 )
      X0 = X(1);
      X = [X0 X];
   elseif ( length(X)>0 )
      X = [X0 X];	
   end


   %Finding Triggers
   if ( length(X) > 0 )
      %Setting Anything < Tresh to zero
      index = find( X < abs(Thresh) );
      X(index) = zeros(1,length(index));

      %Finding Edges in X
      D = diff(X);
      indexD = find( D >= abs(Thresh) );	%Finding Trigger Locations

      %Converting to Trigger Times
      trigger = [trigger M*count+indexD];

      %Finding the Last Element of X. This is placed as the first 
      % element in the new X Array
      X0 = X(length(X));

      %Incrementing Counter
      count = count + 1;
   end
end


%Re-calculating inter-trigger spacing
Dtrigger = diff(trigger);
maxD = max(Dtrigger);

%Finding triple trigger(s)
triple_ind = find(Dtrigger(1:end-1) < 0.25*maxD & Dtrigger(2:end) < 0.25*maxD );
triptrig = trigger(triple_ind);

%Removing additional triggers if extra triple triggers are found
if(isempty(triple_ind))

    disp('Could not find initial triple trigger.');

else %Triple triggers found

   % Multiple triple triggers found. We will assume that the first and 
   % last triple trigger were for the same stimulus. This is the case
   % for the static ripple/moving ripple stimuli, but not for other stimuli
   if( length(triple_ind) == 2 ) 
   %         trigger = trigger( triple_ind(1):(triple_ind(2)-1) );
     trigger = trigger( triple_ind(1):triple_ind(2) );

   %There was only one triple trigger
   else
     trigger = trigger(triple_ind(1):end);
   end

   %Making triple trigger count as one trigger
   trigger = [trigger(1) trigger(4:end)];

end

return;



