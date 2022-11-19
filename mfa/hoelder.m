function hoelder = hoelder(F,int,dim)%(5,6)

 % Generalized Hurst (Hoelder) exponent
 %
 % F is the fluctuation function F(s).
 % int is the detrending interval, for example, int=round(logspace(log10(10),log10(length(x)/4),32)).
 % hoelder = hoelder(F,int,dim) computes the Hurst exponent along the dimension specified by scalar dim.
 
 %hoelder - Hoelder exponent hoelder = hoelder(hoeld,q) of mdfa
 %   hoelder = hoelder(hoeld,q) returns the Scalig exponent hoeld(hoeld(q),q) of
 %   Hoelder exponent hoeld and exponent q via Legandre transform of
 %   hoeld(q) and q.
 %   F is detrended data D
 %   box_size is a vector of detrending range of D.
 %   F is a vector or 2D or 3D matrix. If F is a matrix,
 %   hoeld returns the exponent from each column of the matrix.
 %   See also mdfa,tau,spect
 
 % References
 % 1. Muzy, J., Bacry, E., Arneodo, A.: Multifractal formalism revisited with wavelets. International Journal of Bifurcation and Chaos in Applied Sciences and Engineering 04 (April 1994)
 % 2. Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., Stanley, H.E.: Multifractal detrended ?uctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications 316(1-4) (December 2002) 87?114

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin < 3, dim = 1; end
%F = shiftdim(F);
ndimsx = ndims(F);%(5{)
F = shiftdim(F,dim-1);
sizx = size(F);%(5,6)
F = reshape(F,sizx(1),prod(sizx(2:end)));%(5:6)
if nargin < 2, int = find(F,1,'first'):sizx(1);F = F(int(1):end,:); end
hoelder =-log(F)./log(int(end)/repmat(shiftdim(int(:),1-dim),[sizx(1:dim-1),1,sizx(dim+1:end)]));%(6'')
hoelder = reshape(hoelder,[size(F,1)-1,sizx(2:end)]);%(5)
