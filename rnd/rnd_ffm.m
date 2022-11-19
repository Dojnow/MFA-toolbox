
function y = rnd_ffm(l,g)%(4)

 % rnd_ffm generates correlated time series by modified Fourier filtering method.
 %
 % input parameter l is the length of series; g is the coefficient of correlation.
 %
 %rand_ffm generate random sequence with length L = 2^p, gamma g and type noise_type
 %   y = rand_ffm(L,g) generate fractal noise with length L and gamma g, 0<g<1
 %   y = rand_ffm(L,g,noise_type) generate noise with of type noise_type
 %   if noise_type is 0 or 'w' or 'white' then is generated white noise
 %   if noise_type is 1 or 'f' or 'fractal' then is generated fractal noise
 %   if noise_type is 2 or 'w' or 'brown' then is generated brown noise
 %   See also rand, rand_bmf, rand_ffm, rand_rnd

 % Reference
 % 1. Makse, H.A., Havlin, S., Schwartz, M., Stanley, H.E.: Method for generating long-range correlations for large systems. Phys. Rev. E 53(5) (May 1996) 5445?5449

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

%beta=1-gamma=2*alpha-1=2*H+1;alpha=H+1=1-gamma/2;D=H-2;H=1-gamma/2;0<H<1,'H.Makse'

%if nargin < 3, noise_type = 1; end%(3)
if nargin < 2, g = 1; end%alpha = 0.5;'White_Noise'
%offset = pi/L*2^-3;%0.0001;%(1)
offset = 0.1;%(2,3)
L = 2^l;%(1)
b = (g-1)/2;
%m = -L/2:L/2-1;%('0')%(1,3)
%m = -L/2:L/2;m = m(m~=0);%('~0')%(2.1)
m = [-L/2:-1,1:L/2];%('~0')%(2.2)
q = 2*pi*m/L + offset;%(1:3)
%q0 = find(q==0);q = q(q~=0);%('~0')%(2.2)
%q(q==0) = pi/L*2^-3;%('0+ofs')%(3)
u = randn(1,L);
uq = fft(u);
c = 2*sqrt(pi)/gamma(b+1);
S = c*((q/2).^b).*besselk(b,q);
nq = sqrt(S).*uq;
y = real(ifft(nq));
