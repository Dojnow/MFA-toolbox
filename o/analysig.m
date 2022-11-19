function varargout = analysig(S,dim)%(5.2)

 % Analytic signal otained by Hilbert transform. The amplitude and the phase are derived from an analytical signal.
 %
 % S is a vector or a matrix. If S is a multidimensional array, and  if parameter dim is absemt, then it treats the values along the first non-singleton dimension as vectors. A = mean(S,dim) returns the analytical signal components along the dimension of S specified by scalar dim.
 %
 % Synatax:
 % AP = analysig(S,dim) %Amplitude and phase are concatenatd in one matrix.
 % [A,P] = analysig(S,dim) %A is the amplitude, P is the phase.
 %
 % References
 % 1. D. Gabor, Theory of Communication, Part 1, J. Inst. of Elect. Eng. Part III, Radio and Communication, vol 93, p. 429 1946 (http://genesis.eecg.toronto.edu/gabor1946.pdf)

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
if nargin < 2, dim=1; end%(5.3)
perm = [dim,1:dim-1,dim+1:ndims(S)]; S = permute(S,perm);%(5.2)
sizx = size(S); S = reshape(S,sizx(1),prod(sizx(2:end)));%(5.2)
S = hilbert(S);%(5.1)
S = ipermute(reshape(S,sizx),perm);%(5.2)
Amplitude = abs(S); Phase = angle(S);
switch nargout%(3,4){
    case 1, varargout(1)={[Amplitude,Phase]};%A=[Amplitude,Phase];%(3)
    case 2, varargout(1)={Amplitude};varargout(2)={Phase};%A=Amplitude;B=Phase;%(3)
    %case 3, varargout(1)={[Amplitude,Phase]};varargout(2)={Amplitude};varargout(3)={Phase};
end%(3,4)}
