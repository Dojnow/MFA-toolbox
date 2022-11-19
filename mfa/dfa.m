function [y,int] = dfa(x,int,til)%(7)dfa6(x,int)%(4)

 % Detrended Fluctuation Analysis
 %
 % Y is the fluctuation function F(s). x is a vector or a matrix. int is the detrending interval, for example, int=round(logspace(log10(10),log10(length(x)/4),32)); generates 32 points between 10 and length(x)/4. til is the width of time window.
 %
 % References
 % 1. Peng, C.K., Buldyrev, S.V., Havlin, S., Simons, M., Stanley, H.E., Goldberger, A.L.: Mosaic organization of dna nucleotides. Phys. Rev. E 49(2) (Feb 1994) 1685?1689

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.


x = shiftdim(x);%(5)
[lx1,channels] = size(x);%(2:4)
int = unique(int);%(7)=[int(int(1:end-1)~=int(2:end)),int(end)];%(6)
assert(~(int(1)>int(end) || int(end)>til || til>lx1),'In input arguments')%(7)
%(:6)-'
tinte=til;%(5+)-'
tsegs = floor(lx1/tinte);%(5)-'
% tsegs = lx1-til+1;%(7)="
with1 = tsegs*channels;%(5)-'
x = reshape(x(1:tinte*tsegs,:),tinte,with1);%(5)-'
x = cumsum(x-repmat(mean(x),tinte,1));%(6)-'
x = [x,flipdim(x,1)];%(3)-'
y = zeros(length(int),with1);%(5corr)-'
for s = 1:length(int)%(4)
    nos = floor(tinte/int(s));%(5)-'
    lengthr = int(s)*nos;%(4)-'
    with = nos*with1;%(5)-'
    cs = 1/sqrt(lengthr*2);%(3)-'
    a = [ones(int(s),1),(1:int(s))'];%(4)(all)
    r = reshape(x(1:lengthr,:),int(s),with*2);%(4)-'
    r = r - a*(a\r);
    r = reshape(r,lengthr,with1*2);%(5)-'
    r = [r(:,1:end/2); r(:,1+end/2:end)];%(3sn2,4n2)-'
    for i = 1:with1%(5)-'
        y(s,i) = norm(r(:,i))*cs;%(3n2,4n2)-'
    end
end
y = squeeze(permute(reshape(y,length(int),tsegs,channels),[1,3,2]));%(5)-'