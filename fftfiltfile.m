function fftfiltfile(infile, outfile, h, datatype)
% Filter a large file using the High-speed Overlap-Save method 
% of block convolution using the FFT.
% --------------------------------------------------------------
% fftfiltfile(infile, outfile, h, datatype)
%
% infile = the file of values to be filtered.
%
% outfile = the resulting file after the data in infile 
%    is filtered.
%
% h = filter impulse response. No constraints here, though it
%    should be obvious that a linear phase FIR filter is best.
%
% datatype = type of data in the files. For example, options are
%    'float', 'int16', etc. The default is 'int16'.
%
% fftfiltfile.m optimizes the fft length and the number of blocks
% used in the overlap-save method for best performance.
%
% One way to get the filter h is through Monty's filter design
% functions:
%
% h = lowpass(f2, tw, fs, att, display);
% h = bandpass(f1, f2, tw, fs, att, display);
%
% where f1 and f2 are lower and upper cutoff frequencies,
% tw is the transition width (in Hz), fs is the sampling
% rate of the signal, att is the passband and stopband error,
% and display is 'on' or 'off' to see the filter specs.
%
% Another option is to use the fir1.m or remez.m functions
% for filter design in matlab. These functions also create
% FIR linear phase filters.
%
% This function was created to filter very long files, usually
% over 50 Mb, in a much faster way than implementing the same thing
% using the function filter.m
%
% To eliminate delay accumulation when filtering the long files
% I used the following non-causal steps:
%
% (1) Filter file with fftfiltfile.m
% (2) Flip the file in time.
% (3) Filter the file again with fftfiltfile.m
% (4) Flip the file again in time.
%
% The last output file will be filtered and there will be no delay. 
% 
% caa 2/28/02

if ( nargin < 3 | nargin > 4 )
   error('You need 3 or 4 input args.');
end

if ( nargin == 3 )
   datatype = 'int16';
end


% Find the length of the input file
nx = 0;
fidin = fopen(infile,'r');

while ~feof(fidin)
   temp = fread(fidin, 1024*128, datatype);
   nx = nx + length(temp);
end

fclose(fidin);


% filter specs and form
h = h(:)';
nh = length(h);


%--------- the following is from the function fftfilt.m --------------
fftflops = [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
   149647 319644 680105 1441974 3047619 6422736 13500637 28311786 ...
   59244791];
n = 2.^(1:20); 
validset = find(n>(nh-1));   % must have nfft > (nh-1)
n = n(validset); 
fftflops = fftflops(validset);
% minimize (number of blocks) * (number of flops per fft)
L = n - (nh - 1);
[~,ind] = min( ceil(nx./L) .* fftflops );
nfft = n(ind);
L = L(ind);
%--------------------------------------------------------------------


hk = fft(h,nfft);% Fourier transform of the filter

numblocks = floor((nx+nh-1-1)/(L)); % # of blocks to use

fidin = fopen(infile, 'r');
fidout = fopen(outfile, 'w');

% [nx nh numblocks L nfft]
numblocksamples = (numblocks+1) * L;
numextrasamples = numblocksamples - nx;

for k = 0:numblocks

   if ( k == 0 )
      input = fread(fidin, L, datatype);
      data = [zeros(1,nh-1) input(:)'];
   else
      input1 = data((L+1):nfft);
      input2 = fread(fidin, L, datatype);
      data = [input1(:)' input2(:)'];

      if ( length(data) ~= nfft )
         data = [data zeros(1, nfft-length(data))];
      end

   end

   dk = fft(data);%1. descrete Fourier transform on data segment
   y = real(ifft(dk.*hk));%2. multiply by the filter and then reverse fft
   yout = y(nh:nfft); % get the non-overlapped data
   
   % now output the appropriate number of samples
   if ( k ~= numblocks )
      fwrite(fidout, yout, datatype);
   else
      fwrite(fidout, yout(1:end-numextrasamples), datatype);
   end

end

fclose('all');


% This function was adapted from code intended for signals
% available in the matlab workspace, though not for long
% signals stored in files. The original code, provided in
% the following lines, is from the book by Ingle and Proakis, 
% 'Digital Signal Processing using Matlab.'
%
%
% function [y] =  hsolpsav(x,h,N)
% % High-speed Overlap-Save method of block convolutions using FFT
% % --------------------------------------------------------------
% % [y] = hsolpsav(x,h,N)
% % y = output sequence
% % x = input sequence
% % h = impulse response
% % N = block length (must be a power of two)
% %
% N = 2^(ceil(log10(N)/log10(2));
% Lenx = length(x); M = length(h);
% M1 = M-1; L = N-M1;
% h = fft(h,N);
% %
% x = [zeros(1,M1), x, zeros(1,N-1)];
%
% x = [zeros(1,length(h)-1), x, zeros(1,nfft-1)];
% L = nfft - ( length(h)-1 )
%
% K = floor((Lenx+M1-1)/(L)); % # of blocks
% Y = zeros(K+1,N);
% for k=0:K
%         xk = fft(x(k*L+1:k*L+N));
%         Y(k+1,:) = real(ifft(xk.*h));
% end
% Y = Y(:,M:N)'; y = (Y(:))';


