%'~/My Documents/My Work/work/dfa/m/spect2:4.m','(5)19.02.07,20:03; (6.1)24.04.08,05:31; (6.2)12.05.08,03:10; (6.3)14.05.08,03:21; 17.05.08,08:00,09:53; (6.5)18.12.10,22:39; 09.06.11,21:20;',#2,#3
function [a,f] = spect(tau,q,dim)%(6)

 % spect - Multifractal spetrum f(a) = spect(tau, q) of mdfa
 %
 %   [a, f] = spect(tau, q) returns the Multifractal spetrum f(a) of scaling
 %   exponent tau and exponent q via Legandre transform of tau(hoeld(q),q),
 %   hoeld(q) and q.
 %
 %   tau is a vetor or 2D or 3D matrix. If tau is a matrix,
 %   spect returns the spetrum from each column of the matrix.
 %
 %   See also mdfa, hoelder, tau

 % (\alpha) is the singularity strength or Hoelder exponent, while f(\alpha) denotes the dimension of the subset of the series that is characterized by \alpha.
 % tau is the scaling exponent. q is the index variable, for example, q=[-10:1:10].
 % [a,f] = spect(tau,q,dim) computes the singularity spectrum along the dimension of tau specified by scalar dim.
 %
 % References
 % 1. Muzy, J., Bacry, E., Arneodo, A.: Multifractal formalism revisited with wavelets. International Journal of Bifurcation and Chaos in Applied Sciences and Engineering 04 (April 1994)
 % 2. Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., Stanley, H.E.: Multifractal detrended fluctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications 316(1-4) (December 2002) 87?114
 
 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
 if nargin < 2%(6.1)nargin == 1%(5)
    q = [-2.^(4:-1:-3),0,2.^(-3:4)];
end
if nargin < 3%(6.)
    %dim = 1;%(6.1)
    dim = find(length(q)==size(tau),1,'first');%(6.4)
    if isempty(dim)%(6.4)
        disp('Error'),return
    end
end
%tau = shiftdim(tau);%(5)
ndimsx = ndims(tau);%(6.2)
tau = shiftdim(tau,dim-1);%(6.2)
sizx = size(tau);%(5,6)
tau = reshape(tau,sizx(1),prod(sizx(2:end)));%(5:6.3)
%q = repmat(q',[1,sizx(2:end)]);%(4)
q = repmat(q(:),1,prod(sizx(2:end)));%(5:6.3)
%q = repmat(q',[sizx(1:dim-1),1,sizx(dim+1:end)]);%(6.1err)
%q = repmat(shiftdim(q,find(size(q)~=1)-dim),[sizx(1:dim-1),1,sizx(dim+1:end)]);%(6.1)
%q = repmat(shiftdim(q(:),1-dim),[sizx(1:dim-1),1,sizx(dim+1:end)]);%(6.2)
a = diff(tau)./diff(q);
% a = (tau(3:end)-tau(1:end-2))./(q(3:end)-q(1:end-2));%(6.5.3)
%f = q(2:end,:,:).*a - tau(2:end,:,:);%(4)
f = q(2:end,:).*a - tau(2:end,:);%(5)
%f = q(1:end-1,:).*a - tau(1:end-1,:);%(6.5.2)
% f = q(2:end-1,:).*a - tau(2:end-1,:);%(6.5.3)
a = reshape(a,[sizx(1)-1,sizx(2:end)]);%(6.4)
%a = reshape(a,[sizx(1)-2,sizx(2:end)]);%(6.5)
%a = shiftdim(a,ndimsx-dim);%(6.2)err
a = shiftdim(a,ndimsx-dim+1);%(6.4)
%a = permute(a,[1:dim-1,dim+1:ndimsx,dim]);%(6.3)err
f = reshape(f,[sizx(1)-1,sizx(2:end)]);%(5)
%f = reshape(f,[sizx(1)-2,sizx(2:end)]);%(6.5)
%f = shiftdim(f,ndimsx-dim);%(6.2)err
f = shiftdim(f,ndimsx-dim+1);%(6.4)
%f = permute(f,[1:dim-1,dim+1:ndimsx,dim]);%(6.3)err
