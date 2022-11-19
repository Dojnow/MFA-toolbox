function y = set_cantor(n,p)%(1.3,4)

 % set_cantor generates arbitrary Cantor set.
 %
 % n is the length of set, dependent from pattern. p is the pattern of the set.
 %
 % Example:
 % p=[0.5 0 0.5]; %Monofractal Cantor set. p=[0.1 0 0.6 0 0.3]; %Multifractal Cantor set.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

%function y = cantor(n)%(1:3) n=3^16,(1)21.4s,(2)5.45s,(3)4.95s
%if n<4, y=[1 0 1]; else n=round(n/3);y=[cantor(n),zeros(1,n),cantor(n)]; end%(1.1)
%if n<2, y=1; else n=round(n/3);y=[cantor(n),zeros(1,n),cantor(n)]; end%(1.2)%n=3^n;
%if n<lp+1, y=p; else n=round(n/lp);y=[cantor(n)*p(1),cantor(n)*p(2),...,cantor(n)*p(lp)]; end%(1.3)
%y=1;%(2.2) y=[1 0 1];%(2.1) %length(y)=3^n;%(2)
%y=zeros(1,3^n);y(1)=1;%=[1,zeros(1,3^(n-1))];%(3)
lp=length(p); y=zeros(1,lp^n); p=p(:)'; y(1:lp)=p;%p=shiftdim(p)';%(4)
for i = 1:n-1%(4)=1:n%(2:4)
    %y=[y,zeros(1,3^(i-1)),y];%(2.2)=[y,zeros(1,length(y)),y];%(2.1)
    %y(2*3^(i-1)+1:3^i) = y(1:3^(i-1));%(3)
    lpi=lp^i; lpi1=lpi*lp;%(4)
    %y(1:lpi1)=repmat(y(1:lpi),1,lp).*repmat(p,1,lpi);%(4.1)err
    y(1:lpi1)=reshape(repmat(y(1:lpi),lp,1),1,lpi1).*repmat(p,1,lpi);%(4.2)
    %y(1:lpi1)=repmat(y(1:lpi),1,lp).*reshape(repmat(p,lpi,1),1,lpi1);%(4.3)
end%(2:4)%devil_staircase = cumsum(y);%(1:2)
