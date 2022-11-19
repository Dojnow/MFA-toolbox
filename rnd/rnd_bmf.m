
function y = rnd_bmf(p,a)

 % rnd_bmf generates multifractal time series by binomial multifractal model.
 %
 % input parameter p is the lenght of series; a is the coefficient of multifractality.
 %
 % Reference
 % 1. Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., Stanley, H.E.: Multifractal detrended ?uctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications 316(1-4) (December 2002) 87?114

%random binomial multifractal; a is in(0.5..1)
%rand_bmf generate random sequence with length '2^p' and multifractality 'a'.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

N = 2^p; y = zeros(1,N);
for k = 1:N
    n1 = nnz(de2bi(k));
    y(k) = a^n1*(1-a)^(p-n1);
end


function b = de2bi(varargin)
%DE2BI Convert decimal numbers to binary numbers.
%   B = DE2BI(D) converts a nonnegative integer decimal vector D to a
%   binary matrix B. Each row of the binary matrix B corresponds to one
%   element of D. The default orientation of the binary output is
%   Right-MSB; the first element in B represents the lowest bit.
%
%   In addition to the vector input, three optional parameters can be
%   given:
%
%   B = DE2BI(...,N) uses N to define how many digits (columns) are output.
%
%   B = DE2BI(...,N,P) uses P to define which base to convert the decimal
%   elements to.
%
%   B = DE2BI(...,MSBFLAG) uses MSBFLAG to determine the output
%   orientation.  MSBFLAG has two possible values, 'right-msb' and
%   'left-msb'.  Giving a 'right-msb' MSBFLAG does not change the
%   function's default behavior.  Giving a 'left-msb' MSBFLAG flips the
%   output orientation to display the MSB to the left.
%
%   Examples:
%   >> D = [12; 5];
%
%   >> B = de2bi(D)                 >> B = de2bi(D,5)
%   B =                             B =
%        0     0     1     1             0     0     1     1     0
%        1     0     1     0             1     0     1     0     0
%
%   >> T = de2bi(D,[],3)            >> B = de2bi(D,5,'left-msb')
%   T =                             B =
%        0     1     1                   0     1     1     0     0
%        2     1     0                   0     0     1     0     1
%
%   See also BI2DE.
%   Copyright 1996-2005 The MathWorks, Inc.
%   $Revision: 1.19.4.2 $  $Date: 2004/11/17 20:35:11 $
% Typical error checking.
error(nargchk(1,4,nargin));
% --- Placeholder for the signature string.
sigStr = '';
msbFlag = '';
p = [];
n = [];
% --- Identify string and numeric arguments
for i=1:nargin
   if(i>1)
      sigStr(size(sigStr,2)+1) = '/';
   end;
   % --- Assign the string and numeric flags
   if(ischar(varargin{i}))
      sigStr(size(sigStr,2)+1) = 's';
   elseif(isnumeric(varargin{i}))
      sigStr(size(sigStr,2)+1) = 'n';
   else
      error('Only string and numeric arguments are accepted.');
   end;
end;
% --- Identify parameter signitures and assign values to variables
switch sigStr
   % --- de2bi(d)
   case 'n'
      d		= varargin{1};
	% --- de2bi(d, n)
	case 'n/n'
      d		= varargin{1};
      n		= varargin{2};
	% --- de2bi(d, msbFlag)
	case 'n/s'
      d		= varargin{1};
      msbFlag	= varargin{2};
	% --- de2bi(d, n, msbFlag)
	case 'n/n/s'
      d		= varargin{1};
      n		= varargin{2};
      msbFlag	= varargin{3};
	% --- de2bi(d, msbFlag, n)
	case 'n/s/n'
      d		= varargin{1};
      msbFlag	= varargin{2};
      n		= varargin{3};
	% --- de2bi(d, n, p)
	case 'n/n/n'
      d		= varargin{1};
      n		= varargin{2};
      p  	= varargin{3};
	% --- de2bi(d, n, p, msbFlag)
	case 'n/n/n/s'
      d		= varargin{1};
      n		= varargin{2};
      p  	= varargin{3};
      msbFlag	= varargin{4};
	% --- de2bi(d, n, msbFlag, p)
	case 'n/n/s/n'
      d		= varargin{1};
      n		= varargin{2};
      msbFlag	= varargin{3};
      p  	= varargin{4};
	% --- de2bi(d, msbFlag, n, p)
	case 'n/s/n/n'
      d		= varargin{1};
      msbFlag	= varargin{2};
      n		= varargin{3};
      p  	= varargin{4};
   % --- If the parameter list does not match one of these signatures.
   otherwise
      error('Syntax error.');
end;
if isempty(d)
   error('Required parameter empty.');
end
inType = class(d);
d = double(d(:));
len_d = length(d);
if max(max(d < 0)) || max(max(~isfinite(d))) || (~isreal(d)) || (max(max(floor(d) ~= d)))
   error('Input must contain only finite real positive integers.');
end
% Assign the base to convert to.
if isempty(p)
    p = 2;
elseif max(size(p) ~= 1)
   error('Destination base must be scalar.');
elseif (~isfinite(p)) || (~isreal(p)) || (floor(p) ~= p)
   error('Destination base must be a finite real integer.');
elseif p < 2
   error('Cannot convert to a base of less than two.');
end;
% Determine minimum length required.
tmp = max(d);
if tmp ~= 0 				% Want base-p log of tmp.
   ntmp = floor( log(tmp) / log(p) ) + 1;
else 							% Since you can't take log(0).
   ntmp = 1;
end
% This takes care of any round off error that occurs for really big inputs.
if ~( (p^ntmp) > tmp )
   ntmp = ntmp + 1;
end
% Assign number of columns in output matrix.
if isempty(n)
   n = ntmp;
elseif max(size(n) ~= 1)
   error('Specified number of columns must be scalar.');
elseif (~isfinite(n)) || (~isreal(n)) || (floor(n) ~= n)
   error('Specified number of columns must be a finite real integer.');
elseif n < ntmp
   error('Specified number of columns in output matrix is too small.');
end
% Check if the string msbFlag is valid.
if isempty(msbFlag)
   msbFlag = 'right-msb';
elseif ~(strcmp(msbFlag, 'right-msb') || strcmp(msbFlag, 'left-msb'))
   error('Invalid string msbFlag.');
end
% Initial value.
b = zeros(len_d, n);
% Perform conversion.
%Vectorized conversion for P=2 case
if(p==2)
    [f,e]=log2(max(d)); % How many digits do we need to represent the numbers?
    b=rem(floor(d*pow2(1-max(n,e):0)),p);
    if strcmp(msbFlag, 'right-msb')
        b = fliplr(b);
    end;
else
    for i = 1 : len_d                   % Cycle through each element of the input vector/matrix.
        j = 1;
        tmp = d(i);
        while (j <= n) && (tmp > 0)     % Cycle through each digit.
            b(i, j) = rem(tmp, p);      % Determine current digit.
            tmp = floor(tmp/p);
            j = j + 1;
        end;
    end;
    % If a msbFlag is specified to flip the output such that the MSB is to the left.
    if strcmp(msbFlag, 'left-msb')
        b2 = b;
        b = b2(:,n:-1:1);
    end;
end;
b = feval(inType, b);   % data type conversion
% [EOF]
