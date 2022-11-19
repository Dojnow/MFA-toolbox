function varargout = wtmm_c(x,q,int,til,prec,csel)%(4)

 % wtmm_c is a wrapper function to a C function multifractal at https://archive.physionet.org/physiotools/multifractal/ which implements Wavelet Transform Modilus Maxima method
 %
 % x is a vector or a 2D matrix. q is the index variable, for example, q=[-10:0.2:10].
 % The other input parameters are dummy for compatiblity format with mdfa function.
 % The order of the used Gaussian derivative wavelet is 1.
 %
 % Syntax:
 % Y = wtmm(x,q);
 % [Y,int] = wtmm(x,q);
 % [Y,int,q] = wtmm(x,q);
 %
 % References:
 % 1. https://archive.physionet.org/physiotools/multifractal/
 % 2. Muzy, J., Bacry, E., Arneodo, A.: Multifractal formalism revisited with wavelets. International Journal of Bifurcation and Chaos in Applied Sciences and Engineering 04 (April 1994)

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

p=mfilename('fullpath'); zi=[p,'/zi']; zo=[p,'/zo']; zp='./multifractal'; zp=['"',p,zp(2:end),'"'];%p='oem/wav/wtmm/wtmm';%(20151001221007)
lx=length(x); stp=0.2; q=q(1):stp:q(end); if nargin<4, til=false; end
dlmwrite(zi,[(1:lx)',x(:,1)],' ')%(2)
unix([zp,' "',zi,'" ',num2str(lx),' ',num2str(q(1)),' ',num2str(q(end)+stp),' 1 1 >"',zo,'"']);%(2)
Z = dlmread(zo,' ', 1, 0); Y=zeros(size(Z,1),length(q),size(x,2));%(2)
%int=10.^Z(:,1); Y(:,:,1)=repmat(int,1,length(q)).*10.^Z(:,2:end-1)./repmat(q,length(int),1);%(2)
int=Z(:,1); Y(:,:,1)=Z(:,2:end-1); cond = til==true||isnumeric(til)&&til<40; %(3)'" %(6)'''
if cond, int=10.^int; Y(:,:,1)=repmat(int,1,length(q)).*10.^Y(:,:,1)./repmat(q,length(int),1); end%(6) if til,%(3) if til==true||til<40,%(5)
for i=2:size(x,2)%(2)%1:size(x,2)%(1)
    dlmwrite(zi,[(1:lx)',x(:,i)],' ')
    unix([zp,' "',zi,'" ',num2str(lx),' ',num2str(q(1)),' ',num2str(q(end)+stp),' 1 1 >"',zo,'"']);
    Z = dlmread(zo,' ', 1, 0);
    %Y(:,:,i)=repmat(int,1,length(q)).*10.^Z(:,2:end-1)./repmat(q,length(int),1);%(2) int(:,i)=Z(:,1);%(1)
    Y(:,:,i)=Z(:,2:end-1); if cond, Y(:,:,i)=repmat(int,1,length(q)).*10.^Y(:,:,i)./repmat(q,length(int),1); end %(3)' %(6)" if til,%(3)
end
switch nargout%(4){
    case 1
        varargout(1)=Y;%varargout(1)={[repmat(int,size(Y(1,1,:))),Y]};
    case 2
        varargout(1)={Y}; varargout(2)={int};
    case 3
        varargout(1)={Y}; varargout(2)={int}; varargout(3)={q};
end%(4)}
% delete(zi), delete(zo)%(3)