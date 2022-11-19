function [a1,a0] = fit(X,fitint,l,dim)%(8)

 % Fitting with liner function
 %
 %a1 and a0 are coefficients of liner regression of X, y=a1*x'+a0. fitint is the fitting interval.
 %[a1,a0] = fit(X,fitint,l,dim) fits X along the dimension specified by scalar dim. l is 'lin' - linear or 'log' - logarithmic fitting.

 %   X is a vector or 2D or 3D matrix. If X is a matrix,
 %   fit_l returns the angle from each column of the matrix.
 %   interval is a vector of range of x.
 %   Int is a starting range
 %   See also dfa,mdfa

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin < 4%(7)
    dim = find(length(fitint)==size(X),1,'first');%(7)%dim=1;%(6.1)
    if isempty(dim), disp('Error'),return, end%(7)
end
ndimsx = ndims(X);%(7)
X = shiftdim(X,dim-1);%(7)
sizx = size(X);%(2:6)
%(find(X,1,'first')-1)*prod(sizx(2:end))+1
X = reshape(X,sizx(1),prod(sizx(2:end)));%(5:7)
fl = fitint(isletter(fitint));%(9) fl=fitint(1:3);%(7) fl=fitint(1);%(6)
%if nargin == 1 || nargin == 2 && fl == 'l'%(6)%if nargin == 1%(5)
if nargin == 1 || nargin == 2 && (strcmp(fl,'log') || strcmp(fl,'lin'))%(7)
    fitint = find(X,1,'first'):sizx(1); X = X(fitint(1):end,:);%(5:7)
end%(6,7)
if  nargin == 2 && strcmp(fl,'log') || nargin >= 3 && strcmp(l,'log')%(7)
    X = log(X); fitint = log(fitint);%(5:7)
end%(6,7)
A = [ones(length(fitint),1),fitint(:)]\X;%(5:7)
a0 = reshape(A(1,:),[sizx(2:end),1]); a0 = shiftdim(a0,ndimsx-dim);%(8.1)
a1 = reshape(A(2,:),[sizx(2:end),1]); a1 = shiftdim(a1,ndimsx-dim);%(8.1)
