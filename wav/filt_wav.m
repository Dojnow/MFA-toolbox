function Y = filt_wav(X,wavelet,level,dim)%(6)

 % filt_wav - lowpass filtering with multilevel wavelet decomposition
 % Perform a single-level StationaryWaveletDecomposition.&Reconstruction.
 % D is a vector; If D is a matrix, the filter operates on the columns of D. If D is a multidimensional array, the filter operates on the first nonsingleton dimension or on dimension dim.
 % wavelet is a wavelet name;
 % level - wavelet level decomposition.
%filt_wav - Wavelet fiting of vector x
%   x is a vector
%   lev - level of fiting
%   wav - name of wavlet

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin < 4, dim = 1; end%(6,7)
ndimsx = ndims(X);%(6)
X = shiftdim(X,dim-1);%(6)
sizx = size(X);%(6)
X = reshape(X,sizx(1),prod(sizx(2:end)));%(6)
Y = zeros(sizx(1),prod(sizx(2:end)));%(6)
if nargin < 2, wavelet = 'sym4'; end%(7)
if nargin < 3%(7)
    level = str2double(wavelet(end));%(7)
    if isnan(level), level = 1; end%(6,7)
end%(7)
limit = round(linspace(1,sizx(1),length(level)+1));%(7.3)
for j = 1:length(level)%(7)
    for i = limit(j):limit(j+1);%(7.2)
        [c,l] = wavedec(X(i,:),level(j),wavelet);%(7.3)
        Y(i,:) = wrcoef('a',c,l,wavelet,level(j));%(7.3)
    end%(2j,4level,5i,6i)
end%(1i,2i,4column,5j,6j)
Y = shiftdim(reshape(Y,sizx),ndimsx-dim+1);%(6,7)
