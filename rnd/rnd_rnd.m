function y = rnd_rnd(x,param,noise_type)%(3)
%function y = rnd_rnd(x,noise_type)%(1:2)'17.02.07,22:08'

 % rnd_rnd generates random sequence with length x and type noise_type (param);
 %
 % y = y(x) generate fractal noise with length x
 % y = y(x,noise_type) generate noise with of type noise_type
 % if noise_type is 0 or 'w' or 'white' then white noise is generated
 % if noise_type is 1 or 'f' or 'fractal' then fractal (pink) noise is generated
 % if noise_type is 2 or 'b' or 'brown' then brown noise is generated
 % if noise_type is not specified, then pink noise is generated.
 %
 % Reference
 % 1. Jeong, J., Joung, M.K., Kim, S.Y.: Quantification of emotion by nonlinear analysis of the chaotic dynamics of electroencephalograms during perception of 1/f music. Biological Cybernetics 78 (1998) 217?225 10.1007/s004220050428.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin < 3, noise_type  = param; end%(3.1)
if nargin < 2, noise_type  = 1; end%(1:2)
r = zeros(7,1); y = zeros(1,x);%(3.2)
switch noise_type
    case {0,'w','white'}
        y = randint(x,1,[1,120])';%(2)%y=fix(rand(x)*120)+1;%(1)
    case {1,'f','fractal'}
        %r=fix(rand(20,1)*6)+1;%(1){
        %for i=1:x
        %r([fix(rand(7,1)*20)+1])=fix(rand(7,1)*6)+1;y(i)=sum(r);
        %end%(1)}
        for i = 1:x
            q = 1; %Removes equal rand. numbers in p.
            while q
                p = randint(7,1,[1,20]);
                q = length(p)~=length(unique(p));
            end
            r(p) = randint(7,1,[1,6]); y(i) = sum(r);
        end
    case {2,'b','brown'}
        %y = cumsum([fix(rand*120)+1; fix(rand(x,1)*3)-1]);%(1)
        %y = [0; cumsum(randint(x-1,1,[-1,1]))]';%(2.1)
        y = cumsum([randint(1,1,[1,120]);randint(x-1,1,[-1,1])])';%(2.2)
    otherwise error('In input arguments')%warning('MATLAB:rnd_rnd:Invalid noise_type');
end


function out = randint(varargin)
%RANDINT Generate matrix of uniformly distributed random integers.
%   OUT = RANDINT generates a "0" or "1" with equal probability.
%
%   OUT = RANDINT(M) generates an M-by-M matrix of random binary numbers.
%   "0" and "1" occur with equal probability.
%
%   OUT = RANDINT(M,N) generates an M-by-N matrix of random binary numbers.
%   "0" and "1" occur with equal probability.
%
%   OUT = RANDINT(M,N,IRANGE) generates an M-by-N matrix of random integers.
%
%   IRANGE can be either a scalar or a two-element vector:
%   Scalar : If IRANGE is a positive integer, then the output integer
%            range is [0, IRANGE-1].  If IRANGE is a negative integer,
%            then the output integer range is [IRANGE+1, 0].
%   Vector : If IRANGE is a two-element vector, then the output
%            integer range is [IRANGE(1), IRANGE(2)].
%
%   OUT = RANDINT(M,N,IRANGE,STATE) resets the state of RAND to STATE.
%
%   Examples:
%   ? out = randint(2,3)               ? out = randint(2,3,4)
%   out =                              out =
%        0     0     1                      1     0     3
%        1     0     1                      2     3     1
%
%   ? out = randint(2,3,-4)            ? out = randint(2,3,[-2 2])
%   out =                              out =
%       -3    -1    -2                     -1     0    -2
%       -2     0     0                      1     2     1
%
%   See also RAND, RANDSRC, RANDERR.
%   Copyright 1996-2005 The MathWorks, Inc.
%   $Revision: 1.19.4.1 $  $Date: 2004/11/17 20:35:43 $
% Basic function setup.
error(nargchk(0,4,nargin));
% --- Placeholder for the signature string.
sigStr = '';
m = [];
n = [];
range = [];
state = [];
% --- Identify string and numeric arguments
for i=1:nargin
   if(i>1)
      sigStr(size(sigStr,2)+1) = '/';
   end;
   % --- Assign the string and numeric flags
   if(isnumeric(varargin{i}))
      sigStr(size(sigStr,2)+1) = 'n';
   else
      error('Only numeric arguments are accepted.');
   end;
end;
% --- Identify parameter signitures and assign values to variables
switch sigStr
   % --- randint
   case ''
   % --- randint(m)
   case 'n'
      m		= varargin{1};
	% --- randint(m, n)
	case 'n/n'
      m		= varargin{1};
      n		= varargin{2};
	% --- randint(m, n, range)
	case 'n/n/n'
      m		= varargin{1};
      n  	= varargin{2};
      range = varargin{3};
	% --- randint(m, n, range, state)
	case 'n/n/n/n'
      m		= varargin{1};
      n		= varargin{2};
      range = varargin{3};
      state = varargin{4};
   % --- If the parameter list does not match one of these signatures.
   otherwise
      error('Syntax error.');
end;
if isempty(m)
   m = 1;
end
if isempty(n)
   n = m;
end
if isempty(range)
   range = [0, 1];
end
len_range = size(range,1) * size(range,2);
% Typical error-checking.
if (~isfinite(m)) || (~isfinite(n))
   error('Matrix dimensions must be finite.');
elseif (floor(m) ~= m) || (floor(n) ~= n) || (~isreal(m)) || (~isreal(n))
   error('Matrix dimensions must be real integers.');
elseif (m < 0) || (n < 0)
   error('Matrix dimensions must be positive.');
elseif (length(m) > 1) || (length(n) > 1)
   error('Matrix dimensions must be scalars.');
elseif len_range > 2
   error('The IRANGE parameter should contain no more than two elements.');
%elseif max(max(floor(range) ~= range)) | (~isreal(range)) | (~isfinite(range))%(0)
elseif any(floor(range) ~= range) || ~isreal(range) || ~any(isfinite(range))%(APL)
   error('The IRANGE parameter must only contain real finite integers.');
end
% If the IRANGE is specified as a scalar.
if len_range < 2
    if range < 0
        range = [range+1, 0];
    elseif range > 0
        range = [0, range-1];
    else
        range = [0, 0];    % Special case of zero range.
    end
end
% Make sure IRANGE is ordered properly.
range = sort(range);
% Calculate the range the distance for the random number generator.
distance = range(2) - range(1);
% Set the initial state if specified.
if ~isempty(state)
   rand('state', state);
end
% Generate the random numbers.
r = floor(rand(m, n) * (distance+1));
% Offset the numbers to the specified value.
out = ones(m,n)*range(1);
out = out + r;
% [EOF] randint.m
