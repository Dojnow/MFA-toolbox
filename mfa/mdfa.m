function [Y,int,q] = mdfa(x,q,int,til,prec,varargin)%(7.4)

% Multifractal Detrended Fluctuation Analysis implements modified algorithm with sliding windows (segments)
 %
 % Y is the fluctuation function Fq(s,t).
 % x is a vector or a 2D-matrix.
 % q is the index variable can take any real variable, for example, q=[-10:1:10].
 % int the is detrending interval, for example, int=round(logspace(log10(10),log10(length(x)/4),32)).
 % til is the width of time sliding window. If til = int then the entire time series is retrieved.
 % prec is the precision, single or double.
 %
 % References
 % 1. Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., Stanley, H.E.: Multifractal detrended ?uctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications 316(1-4) (December 2002) 87–114
 % 2. Podobnik, B., Stanley, H.E.: Detrended cross-correlation analysis: A new method for analyzing two nonstationary time series. Phys. Rev. Lett. 100 (Feb 2008) 084102

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

x = shiftdim(x);%(6.2)(0;4)(1;4.2)(2;4.2)(3:6;)
[lx1,channels] = size(x);
if nargin<6, cores=feature('numcores'); else cores=varargin{1}; end%(7.9.2)
if nargin<5, prec = 'double'; end, x = cast(x,prec);%(7.3)
if nargin<4, til = lx1; end%(7.4)
assert(~(int(1)>int(end) || int(end)>til || til>lx1),'In input arguments')%(6.9)
tsegs = lx1-til+1;%(6)
int = unique(int);%(7.2?)=[int(int(1:end-1)~=int(2:end)),int(end)];%(0;6)(1;6)(2;6)(3;8)(4;7)
q0 = q==0;%q0 = find(q==0);%(0;3.4)(1;4)(2;4)(3;2:7)(4;2:5)
qix = find(q~=0);%(0;4)(1;5)(2;6)(3;7)(4;5)
q = shiftdim(q);%(1;5)(2;6)(3;7)(4;5)
q_1 = 1./q;%(6) =-1./q;%(3;7)(4;5)(5;)
Y = zeros(length(int),length(q),channels,tsegs,prec);%(7;) @mdfa %
x = cumsum(x-repmat(mean(x),lx1,1));%(4;K,O,D,P)(6;)
preci = prec;%(7) preci='single';%(7.5)
R = zeros(int(end),lx1-int(1)+1,channels,preci);%(7)
if cores==1 %(7.9){for
    for i = 1:channels%(6) =1:tsegs*channels%(D,Dt,P,Pt)%'
        R(:,:,i) = hankel(x(1:int(end),i),x([int(end):lx1,ones(1,int(end)-int(1))],i));%(D,Dt,P,Pt)%' %
    end%' %(6.2)Rout}
    if ~(til<lx1)
        for s = 1:length(int)%(4)
            nos = lx1-int(s)+1;%(6) =tinte+1-int(s);%(D,Dt,P,Pt,4,5)
            a = [ones(int(s),1),(1:int(s))'];%(4)
            with = nos*channels;%(6.2)
            sqrt_s = 1/sqrt(int(s)+1);%(int(s)-1)^-0.5;(6.2)
            q0_ = 1/(til-int(s)+1);%(6) = 1/(nos);%(P,Pt)
            q_ = q0_.^q_1;%(6) = nos.^q_1;%(5)
            r = reshape(R(1:int(s),1:nos,:),int(s),with);%(6.2){Rout} %
            r = r - a*(a\r);
            for i = 1:with%withi%(6.4)
                r(i) = norm(r(:,i))*sqrt_s;%(6.4.2)
            end%(6.4)
            r = reshape(r(1:with),nos,channels);%(6.4.2)
            for j = 1:channels
                for t = 1:tsegs
                    f = r(t:t+til-int(s),j);%(6.4.2)
                    Y(s,q0,j,t) = exp(q0_*sum(log(f)));%(Dt,Pt)
                    for i = qix%(2:4)
                        Y(s,i,j,t) = q_(i)*norm(f,q(i));%(Pt,6.1)
                    end
                end
            end%for j
        end
    elseif any(q0)
        q_1=repmat(q_1(qix),1,channels); q_2=q'; q_2(q0)=1;%(7.6.7)
        for s = 1:length(int)%(4)
            nos = lx1-int(s)+1;%(6) =tinte+1-int(s);%(D,Dt,P,Pt,4,5)
            a = [ones(int(s),1),(1:int(s))'];%(4)
            with = nos*channels;%(6.2)
            sqrt_s = 1/sqrt(int(s)+1);%(int(s)-1)^-0.5;(6.2)
            q0_ = 1/(til-int(s)+1);%(6) = 1/(nos);%(P,Pt)
            segs = til-int(s)+1;%(7.6)
            r = reshape(R(1:int(s),1:nos,:),int(s),with);%(6.2){Rout} %
            r = r - a*(a\r);
            fq = zeros(with,length(q));%(7.6.2:7)
            for i = 1:with%(7.6.2:7)
                fq(i,:) = (norm(r(:,i))*sqrt_s).^q_2;%(7.6.2:3)
            end%(7.6.2:7)
            fq(:,q0)=log(fq(:,q0));%(7.6.2:7) %{any(q0)}
            fq = permute(reshape(fq,nos,channels,length(q)),[3 2 1]);%(7.6.7)
            fqsum = sum(fq(:,:,1:segs-1),3) + fq(:,:,1);%(7.6.5:7)
            for t = 1:tsegs
                fqsum = fqsum - fq(:,:,t) + fq(:,:,t+segs-1);%(7.6.6:7)
                Y(s,q0,:,t) = exp(q0_*fqsum(q0,:));%(7.6.7) %{any(q0)}
                Y(s,qix,:,t) = (q0_*fqsum(qix,:)).^q_1;%(7.6.7) %{any(q0)}
            end
        end
    else
        q_1=repmat(q_1(qix),1,channels); q_2=q'; q_2(q0)=1;%(7.6.7)
        for s = 1:length(int),disp(s)%(4)
            nos = lx1-int(s)+1;%(6) =tinte+1-int(s);%(D,Dt,P,Pt,4,5)
            a = [ones(int(s),1),(1:int(s))'];%(4)
            with = nos*channels;%(6.2)
            sqrt_s = 1/sqrt(int(s)+1);%(int(s)-1)^-0.5;(6.2)
            q0_ = 1/(til-int(s)+1);%(6) = 1/(nos);%(P,Pt)
            segs = til-int(s)+1;%(7.6)
            r = reshape(R(1:int(s),1:nos,:),int(s),with);%(6.2){Rout} %
            r = r - a*(a\r);
            fq = zeros(with,length(q));%(7.6.2:7)
            for i = 1:with%(7.6.2:7)
                fq(i,:) = (norm(r(:,i))*sqrt_s).^q_2;%(7.6.2:3)
            end%(7.6.2:7)
            fq = permute(reshape(fq,nos,channels,length(q)),[3 2 1]);%(7.6.7)
            fqsum = sum(fq(:,:,1:segs-1),3) + fq(:,:,1);%(7.6.5:7)
            for t = 1:tsegs
                fqsum = fqsum - fq(:,:,t) + fq(:,:,t+segs-1);%(7.6.6:7) {for}
                Y(s,:,:,t) = (q0_*fqsum).^q_1;%(7.7) {~any(q0),for}
            end
        end
    end
else %(7.9){parfor}
    %matlabpool(cores)%(7.8) {parfor}
    parfor i = 1:channels%(6) =1:tsegs*channels%(D,Dt,P,Pt)%'
        R(:,:,i) = hankel(x(1:int(end),i),x([int(end):lx1,ones(1,int(end)-int(1))],i));%(D,Dt,P,Pt)%' %
    end%' %(6.2)Rout}
    if ~(til<lx1)
        for s = 1:length(int)%(4)
            nos = lx1-int(s)+1;%(6) =tinte+1-int(s);%(D,Dt,P,Pt,4,5)
            a = [ones(int(s),1),(1:int(s))'];%(4)
            with = nos*channels;%(6.2)
            sqrt_s = 1/sqrt(int(s)+1);%(int(s)-1)^-0.5;(6.2)
            q0_ = 1/(til-int(s)+1);%(6) = 1/(nos);%(P,Pt)
            q_ = q0_.^q_1;%(6) = nos.^q_1;%(5)
            r = reshape(R(1:int(s),1:nos,:),int(s),with);%(6.2){Rout} %
            r = r - a*(a\r);
            f2 = zeros(1,with);%(6.4.1) (7.11){parfor}
            parfor i = 1:with%(7.8) {parfor}
                f2(i) = norm(r(:,i))*sqrt_s;%(6.4.1) (7.8){parfor}
            end%(6.4)
            f2 = reshape(f2,nos,channels);%r=(7.8){parfor}
            for j = 1:channels
                for t = 1:tsegs
                    f = f2(t:t+til-int(s),j);%(6.4.1)
                    %Y(s,q0,j,t) = exp(q0_*sum(log(f)));%(Dt,Pt)
                    parfor i = 1:length(q)%(7.11) {any(q0),parfor}
                        Y(s,i,j,t) = q_(i)*norm(f,q(i));%(Pt,6.1)
                    end
                    Y(s,q0,j,t) = exp(q0_*sum(log(f)));%(7.11.3) {any(q0),parfor}
                end
            end%for j
        end
    else
        q_1=repmat(q_1(qix),1,channels); q_2=q'; q_2(q0)=1;%(7.6.7)
        for s = 1:length(int)%(4)
            nos = lx1-int(s)+1;%(6) =tinte+1-int(s);%(D,Dt,P,Pt,4,5)
            a = [ones(int(s),1),(1:int(s))'];%(4)
            with = nos*channels;%(6.2)
            sqrt_s = 1/sqrt(int(s)+1);%(int(s)-1)^-0.5;(6.2)
            q0_ = 1/(til-int(s)+1);%(6) = 1/(nos);%(P,Pt)
            segs = til-int(s)+1;%(7.6)
            r = reshape(R(1:int(s),1:nos,:),int(s),with);%(6.2){Rout} %
            r = r - a*(a\r);
            fq = zeros(with,length(q));%(7.6.2:7)
            parfor i = 1:with%(7.6.2:7)
                fq(i,:) = (norm(r(:,i))*sqrt_s).^q_2;%(7.6.2:3)
            end%(7.6.2:7)
            fq = permute(reshape(fq,nos,channels,length(q)),[3 2 1]);%(7.6.7)
            parfor t = 1:tsegs
                Y(s,:,:,t) = (q0_*sum(fq(:,:,t:segs+t-1),3)).^q_1;%(7.7) {~any(q0),parfor}
            end
        end
    end
    %matlabpool close%(7.8) {parfor}
end
