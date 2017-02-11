function DispSTFT(c, param)
% Display the spectrogram in the specified range
% this file is modified from 'displaySTFT.m' originally by 
% I. W. Selesnick, Polytechnic Inst. of NYU
%
% INPUT
%   c : STFT coefficients (2D array)
%   param.fs : sampling rate
%   param.N : length of FFT
%   param.hop : Hop-size
%   param.Fr : display frequency range (in Hz)
%   param.clim : for 'Clim' of gca
%
% Ilker Bayram, Istanbul Technical University, 2012.

Fd = 1 + round(param.Fr * param.N / param.fs);

cdb = 20 * log10( abs( c(Fd(1):1:Fd(2),:) ) );
imagesc([0 size(c,2) * param.Hop / param.fs], param.F / 1000, cdb);

set(gca,'Clim',clim);

xlabel( 'Time (seconds)' );
ylabel( 'Frequency (kHz)' );