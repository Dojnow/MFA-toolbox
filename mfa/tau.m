%'~/My Documents/My Work/work/dfa/m/tau1:6.m', '13.11.06,16:22; (5) 24.04.08,04:48; (6.1) 12.05.08,02:18; (6.2) 17.05.08,05:19,10:03; (6.3) 03.11.10,04:06;'
function tau = tau(hoeld,q,dim)%(6)
 % Scaling exponent tau
 %
 % tau - tau is the scalar exponent. q is the index variable, for example, q=[-10:1:10].
 % tau = tau(hoeld,q,dim) computes the tau function along the dimension of hoeld specified by scalar dim.
 %
 %   tau = tau(hoeld, q) returns the Scalig exponent tau(hoeld(q),q) of
 %   Hoelder exponent hoeld and exponent q via Legandre transform of
 %   hoeld(q) and q.
 %
 %   hoeld is a vetor or 2D or 3D matrix. If hoeld is a matrix,
 %   tau returns the exponent from each column of the matrix.
 %
 %   See also mdfa, hoelder, spect
 %
 % References
 % 1. Muzy, J., Bacry, E., Arneodo, A.: Multifractal formalism revisited with wavelets. International Journal of Bifurcation and Chaos in Applied Sciences and Engineering 04 (April 1994)
 % 2. Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., Stanley, H.E.: Multifractal detrended ?uctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications 316(1-4) (December 2002) 87?114

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin < 2%(6)nargin == 1%(5)
    q = [-2.^(4:-1:-3),0,2.^(-3:4)];
end
if nargin < 3%(6.)
    %dim = 1;%(6.1)
    dim = find(length(q)==size(hoeld),1,'first');%(6.2)
    if isempty(dim)%(6.2)
        disp('Error'),return
    end
end
%(6.3)%hoeld = shiftdim(hoeld);%(6.2)err
sizx = size(hoeld);
%tau = hoeld.*repmat(q',[1,sizx(2:end)]) - 1;%(5)
%tau = hoeld.*repmat(q(:),[sizx(1:dim-1),1,sizx(dim+1:end)]) - 1;%(6.1)err
tau = hoeld.*repmat(shiftdim(q(:),1-dim),[sizx(1:dim-1),1,sizx(dim+1:end)])-1;%(6.2)
