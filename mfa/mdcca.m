function [Y,int,q,csel] = mdcca(x,q,int,til,prec,csel)%(7.2)

 % Multifractal Detrended Crosscorelation Fluctuation Analysis
 %
 % Y is the fluctuation function Fq(s). X is a vector,  a 2D-matrix or a 3D-matrix.
 % q is the index variable can take any real variable, fot example, q=[-10:1:10].
 % int is the detrending interval, for example, int=round(logspace(log10(10),log10(length(X)/4),32)).
 % til is the width of time sliding window. prec is the precision, single or double.
 % csel selects which time series (channels) of matrix will be crosscorelated.
 % For example, csel=[1 0 1; 0 0 0; 0 1 0] means that only X(:,1,1), X(:,1,3), X(:,3,2) are crosscorelated each other.
 %
 % References
 % 1. Podobnik, B., Stanley, H.E.: Detrended cross-correlation analysis: A new method for analyzing two nonstationary time series. Phys. Rev. Lett. 100 (Feb 2008) 084102

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

x = shiftdim(x);%(6.2)(4;)(3;)(2;4.2)(1;4.2)(0;4)
[lx1,channels] = size(x);
if nargin<6, csel = triu(true(channels),1); end, csel = logical(csel);%(7.2)
if nargin<5, prec = 'double'; end, x = cast(x,prec);%(7.1)
if nargin<4, til = lx1; end%(7.2)
assert(~(int(1)>int(end) || int(end)>til || til>lx1),'In input arguments')%(6.9)
assert(~any(size(csel) > channels),'In input arguments')%(7.6)
%tinte = lx1;tsegs = 1;til = tinte;%(6.)
tsegs = lx1-til+1;%(6)
x = cumsum(x-repmat(mean(x),lx1,1));%(4;K,O,D,P)(6;)
%tinte = til;%(A;)
int = unique(int);%(7.2?)
q0 = q==0;
qix = find(q~=0);%(4;5)
q = shiftdim(q)';%(6.6)
lq = length(q);%(6)
q_1 = 1./q(qix);%(6.7) = 1./q;%(6)
q_2 = q(qix)/2;%(6.7) = q/2;%(6)
%tsegs = floor(lx1/tinte);%(4;D,Dt,P,Pt)(3;6)(2;5)(1;5)(0;4)
% with1 = tsegs*channels;%3;2;(1;5)(0;4)
[ci,cj] = find(csel);%(C;C)
Y = zeros(length(int),length(q),max(ci),max(cj),tsegs,prec);%(7.4)
R = zeros(int(end),lx1-int(1)+1,channels,prec);%(7)
for i = 1:channels%(6) =1:tsegs*channels%(D,Dt,P,Pt)%'
% parfor i = 1:channels%(7.7) {parfor}
	R(:,:,i) = hankel(x(1:int(end),i),x([int(end):lx1,ones(1,int(end)-int(1))],i));%(D,Dt,P,Pt)%' %
end%' %(6.2)Rout}
for s = 1:length(int)%(4)
	nos = lx1-int(s)+1;%(6) =tinte+1-int(s);%(D,Dt,P,Pt,4)
	a = [ones(int(s),1),(1:int(s))'];%(4)
	%with = nos*channels;
	segs = til-int(s)+1;%(6.6)
	ints_1 = 1/(int(s) - 1);%(5)|(int(s)+1)%(7.8+)|1/int(s);%(C4)
	fq = zeros(nos+1,lq);%(6.5,7;7.5.2,.8.2)
	q_ = 1/segs; q0_ = q_*0.5;%(6)
	r = reshape(R(1:int(s),1:nos,:),int(s),nos*channels);%(6.2){Rout} %
% 	r = reshape(R,int(s),nos*channels);%(6.3){Rin},(7){Radi} %
	r = abs(r - a*(a\r));%(7.5) r = r - a*(a\r);
    r = reshape(r,int(s),nos,channels);%(6.4:5)
	for l = 1:length(ci)%(7.6)=lcseli%() =1:length(csel)%(C2,4)
		for t = 1:nos%(6.3:5) 1:til-1%til+1-int(s)%(6.2) =1:til%(6.1) i=nosi%(C4)
			fq(t,:) = r(:,t,ci(l))'*r(:,t,cj(l))*ints_1;%(7.6) =r(:,t,cix)'*r(:,t,cjx)*ints_1;%(7.8){parfor}
		end
		fq(1:nos,q0) = log(fq(1:nos,q0));%(7.6.2)q0 %fq()=%(7.8)~any(q0)
		fq(1:nos,qix) = fq(1:nos,qix).^repmat(q_2,nos,1);%(7.6.2)
		fqsum = sum(fq(1:segs,:));%(6.5:7) f0sum = sum(f0(1:segs));%(6.5,6)
		for t = 1:tsegs%(6) %for% ~t
		%parfor t = 1:tsegs%(7.7) {parfor}
			Y(s,q0,ci(l),cj(l),t) = exp(q0_*fqsum(q0));%(7.4) {for} %Y() ~any(q0)
			Y(s,qix,ci(l),cj(l),t) = (q_*fqsum(qix)).^q_1;%(7.4) {for}
			fqsum = fqsum - fq(t,:) + fq(t+segs,:);%(6.7) {for} %fqsum= ~any(q0)
        end
	end
end