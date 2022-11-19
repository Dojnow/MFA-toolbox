function [PS,f] = psdm(D,nfft,fs)

 % Power Spectral Density estimate (psd) from the rows of a 2D matixes.
 %
 % D is 2D matrix. PS containing the Power Spectral Density of the elements of each column of D. fs is sampling frequency.
 %
 % fs is sampling frequency in Hz
 % nfft is number of FFT points
 
 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

%channels=size(D,2); %channels=17;%(1) %fs=256;
sizx = size(D);%(2)
D = reshape(D,sizx(1),prod(sizx(2:end)));%(2)
[PS(:,1),f] = psd(D(:,1),nfft,fs);%(2)
PS=[PS(:,1),zeros(length(f),prod(sizx(2:end))-1)];%(3) =zeros(,prod(sizx(2:end)))%(2)
for i = 2:prod(sizx(2:end))%(2) =1:size(D,2) %1:min(size(D))%(1)
    PS(:,i) = psd(D(:,i),nfft,fs);
end
%PS = reshape(PS,[length(f),sizx(2:end)]);%(2)
