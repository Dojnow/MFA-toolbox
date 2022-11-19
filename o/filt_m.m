function Y = filt_m(D,fname,fo,ff,fs,ftype)%(5)
%function Y = filt_butter(D,fo,ff,fs)%(4) %Y=filt_butter(D,ff,fs)%(1:3)

 % filtm_m - Multiband filtering;
 %
 % D is a vector; If D is a matrix, filter operates on the columns of D. If D is a multidimensional array, filter operates on the first nonsingleton dimension.
 % fname - filter name: 'besel','butter',... ; fo - filter order;
 % ff - cutoff frequency(ies) (Hz); fs - sampling frequency; ftype - 'low','high','bandpass','stop';
 %
 % Example:
 % rndn=rand(2^10,1); fname='butter'; fo=[4 6 6 10 10]; ff=[1 3; 3 7; 8 12; 18 30; 35 40]; fs=256; ftype='bandpass'; rndnf=filt_m(rndn,fname,fo,ff,fs,ftype); [ps,fs]=psdm(rndnf,[],256); loglog(fs,ps)

 %%filtm_butter - Multiband filtering; Butterworth filter
 % D - nDim matrix; fname - 'besel','butter',... ; fo - filter order;
 % ff - frequency(ies); fs - sampling frequency; ftype - 'low','bandpass','high','stop';
 %fs=256; ff=[4 1 3; 6 3 7; 6 8 12; 10 18 30; 10 35 40];%(1:3)
  %[lx1,channels]=size(D);%(3) %channels=size(D,2); %channels=17;
 
 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin<6 && isscalar(ff), ftype='low'; elseif nargin<6 && isvector(ff), ftype='bandpass'; end%(5)
Y = zeros([size(D),length(fo)]);%(4:5) Y=zeros(lx1,size(ff,1),channels);%(3) %Y=D;%(2)
for i = 1:length(fo)%(4:5) i=1:size(ff,1)%(3) % i=1:5
    %[b,a]=butter(fo(i),ff(i)/(fs/2));%(4) %[b,a]=butter(ff(i,1),[ff(i,2),ff(i,3)]/(fs/2));%(1:3)
    [b,a] = feval(fname,fo(i),ff(i,:)/(fs/2),ftype);%=builtin()%(5)
    %Y=[Y,filtfilt(b,a,D(:,1:channels))];%(2) %D=[D,filtfilt(b,a,D(:,1:channels))];%(1)
    Y(:,:,i) = filtfilt(b,a,D);%(4:5) %Y(:,i,:)=filtfilt(b,a,D);%(3)
end%, Y = squeeze(Y);%(5) %,Y=D;%(1)

