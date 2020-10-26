%
%function [H] = bandpass(f1,f2,TW,Fs,ATT,Disp)
%	
%	FILE NAME       : Band Pass Filter Design
%	DESCRIPTION 	: Impulse response of optimal Band Pass filter 
%                     as designed by Roark/Escabi.
%
%	H               : Calculated impulse response
%	f1              : Lower cutoff frequency (in Hz)
%	f2              : Upper cutoff frequency (in Hz)
%	TW              : Transition width (in Hz)
%	Fs              : Sampling Frequency (in Hz)
%	ATT             : Passband and stopband error (Attenuationin dB)
%	display         : optional -> 'n' to turn display off
%                     default  -> 'y'
%
% (C) Monty A. Escabi, Last Edit August 2006
%
function [H] = bandpass(f1,f2,TW,Fs,ATT,display)

%Preliminaries
if ( nargin==5 )
	display = 'y';
end

%Designing band pass filter 
if f1==0
	H=lowpass(f2,TW,Fs,ATT,'off');
else
	H1=lowpass(f1,TW,Fs,ATT,'off');
	H2=lowpass(f2,TW,Fs,ATT,'off');
	H=H2-H1;
end

%Display
M=2^ceil(log2(length(H)))*8;
if strcmp(display,'y')
	subplot(211)
	HH=abs(fft(H,M));
	plot((1:M)/M*Fs,20*log10(HH),'r'),axis([0 Fs/2 min(20*log10(HH)) 0])
	ylabel('Power Spectrum (dB)'),xlabel('Frequency (Hz)')
	subplot(212)

	plot((1:M)/M*Fs,HH),axis([f1 f2 1-10^(-ATT/20) 1+10^(-ATT/20)])
	ylabel('Linear Passband Amplitude'),xlabel('Frequency (Hz)')
end
