function Y = rnd_perm(X,dim)%(4)

 % rnd_perm - random permutation of elements of matrix X.
 %
 % X is a vector or a matrix. If X is a matrix, rnd_perm(X) returns a row vector containing the permutations (shuffling) of the elements of each column of X.
 % If X is a multidimensional array, rnd_perm(X) is the permutations of the elements along the first nonsingleton dimension of X.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

%function Y = shuffle(X)%(1:3)
if nargin == 1, dim = 1; end%(4)
%X = shiftdim(X);%(4)err
perm = [dim,1:dim-1,dim+1:ndims(X)]; X = permute(X,perm);%(5)
sizx = size(X); X = reshape(X,sizx(1),prod(sizx(2:end))); %(2:5)'%(5)"
%X = [X, X(repmat(randperm(size(X,1))',1,size(X,2)))];%(1)
%X = reshape(X,sizx(1),prod(sizx(2:end)));%(2)
%X = X(repmat(randperm(sizx(1))',1,sizx(2)));%(2)
%Y = reshape(X,[sizx(2:end),2]);%(2)
%Y = X(repmat(randperm(sizx(1))',1,sizx(2:end)));%(3)
%Y = X(repmat(randperm(sizx(dim))',[sizx(1:dim-1),1,sizx(dim+1:end)]));%(4)err
Y = X(randperm(sizx(dim)),:); Y = ipermute(reshape(Y,sizx),perm);%(5)
