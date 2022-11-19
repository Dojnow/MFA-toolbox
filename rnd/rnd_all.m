function y = rnd_all(ver,lx,param,type)

 % Function rnd_all is a wrapper function of all rnd_* functions for generation of (multi)fractal long-range correlation series and surrogate data.
 %
 % ver is a type of generating method: 'bmf' - Binomial multifractal model; 'ffm' - Fourier filtering method; 'perm' - Random permutation of vector; 'rnd' - generated sequence of numbers selected by tossing imaginary dice; 'cantor' - generates cantor set with arbitrary pattern.
 %lx is the lenght of generated time series.
 %param - coefficient of correlation or multifractality
 %type - type of noise: white, pink, brown.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

assert(~(isscalar(lx) && lx<0),'In input arguments')%(2)%assert(~(lx<0),'In...')%(1)
%switch ver
%    case 'bmf', y = rnd_bmf(lx,param);
%    case 'ffm', y = rnd_ffm(lx,param);
%    case 'perm', y = rnd_perm(lx,param);
%    case 'rnd', y = rnd_rnd(lx,param,type);%=rnd_rnd(lx,param);
%    case 'cantor', y = set_cantor(lx,param);
%    otherwise error('In input arguments')%
%end%(1:2)
st=[repmat('rnd',4,1);'set']; y=feval([st(strcmp(ver,{'bmf','ffm','perm','rnd','cantor'}),:),'_',ver],lx,param);%(3)
if ~any(strcmp(ver,{'perm','rnd'}))%(2)%if ~strcmp(ver,'rnd')%(1)
	if nargin < 4, type = 1; end%(2)
    switch type
        case {1,'f','fractal'}
        case {2,'b','brown'}, y = cumsum(y);
        otherwise error('In input arguments')%warning('MATLAB:rnd_ffm:Invalid noise_type');
    end
end
