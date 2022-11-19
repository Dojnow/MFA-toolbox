
function [Y,int,q] = mdfast(x,q,int,til,prec,varargin)%(201309240632.01)

 % Multifractal Detrended Fluctuation Analysis - fast execution and low memory requirements variant of mdfa implements original algorithm with nonoverlaping windows (segments)
 %
 % Y is the fluctuation function Fq(s).
 % x is a vector or a 2D-matrix.
 % q is the index variable can take any real variable, for example, q=[-10:1:10].
 % int is the detrending interval, for example, int=round(logspace(log10(10),log10(length(x)/4),32)).
 % til, prec and varargin are dummy parameters for compatibility with input parameters of the function mdfa.
 % prec the is precision, single or double.
 %
 % References
 % 1. Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., Stanley, H.E.: Multifractal detrended ?uctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications 316(1-4) (December 2002) 87?114

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if nargin<6, csel=1; else csel=varargin{1}; end%(201510011504)
[Y,x,int,q,q0,q_1,qix,channels,csel,tinte,tsegs,with1] = mfa_args(x,q,int,til,csel);%()
x = [x,flipdim(x,1)];%(K,P)
% matlabpool%(8)
%for s = int%(1:3)
for s = 1:length(int)%(4:6)
    %a = [ones(s,1),(1:s)'];%(1:3)
    a=[ones(int(s),1),(1:int(s))'];%(4:6;dfa1)%
    %a=[ones(int(s),1),(1:int(s))',(1:int(s))'.^2];%(dfa2)
    %nos = floor(lx1/s);%(1:3)
    %nos = floor(lx1/int(s));%(4,5)
    nos = floor(tinte/int(s));%(6)
    nos2 = nos*2;%(7) %
    nosi2 = 1:nos2;%(7)
    with = nos*with1;
    %f = zeros(1,nos);%(O)
    f = zeros(1,nos*2);%(K,P)
    %sqrt_s = 1/sqrt(s);%(1:3)
    sqrt_s = 1/sqrt(int(s));%(int(s)-1)%(4,5)
    q0_ = 1/(nos*2);%(O,P)
    %q0_ = 1/(nos*4);%(K)
    %q_ = (nos).^(-1./q);%(O1.2,2,3,P')err
    %q_ = (nos*2).^(-1./q);%(K1.2,2,3,P")cor
    q_ = q0_.^q_1;%(7.2) = (nos*2).^q_1;%(7)
    %r = reshape(x(1:s*nos,:),s,with);%(O)
    %r = reshape(x(1:s*nos,:),s,with*2);%(K,P,1:3)
    r = reshape(x(1:int(s)*nos,:),int(s),with*2);%(K,P,4:6)
    r = r - a*(a\r);
    for j = 1:channels
        for t = 1:tsegs
            %index = (1:nos)+(j-1)*nos;%(O)
            %index = [1:nos,with+(1:nos)]+(j-1)*nos;%(K,P)
            index = [1:nos,(1:nos)+with]+((t-1)+(j-1)*tsegs)*nos;%(6)
            %rix = [r(:,index(i));r(:,index(i)+with)];%(O12)
            for i = nosi2%1:nos2%1:length(index)%(O,K,P)
            %parfor i = nosi2%(8)
                f(i) = norm(r(:,index(i)))*sqrt_s;%(K1)
                %f(i) = sqrt(r(:,index(i))'*r(:,index(i))/s);%(K2)
                %f(i) = norm(rix)/sqrt(s*2);%(O1)
                %f(i) = sqrt((rix'*rix)/(s*2));%(O2)
            end
            %for i = 1:length(q)%(1)
            %Y(s,i,j) = (mean(f.^qh(i)))^(1/q(i));%(.1)
            %Y(s,i,j) = q_(i)*sqrt(norm(f,q(i)/2));%(.2)
            %end
            %if q0%(2)
            %ix = ix0;%(.1)
            %Y(s,q0,j) = exp(q0_*sum(log(f)));%(.1,.2)
            %else%(.1)
            %ix = ix1;
            %end
            %Y(s,q0,j) = exp(q0_*sum(log(f)));%(3:5)
            for i = qix%(2:6) parfor i=qix(1):qix(end)%(8)mc
                %Y(s,i,j) = q_(i)*sqrt(norm(f,q(i)/2));%(O,K)
                %Y(s,i,j) = q_(i)*norm(f,q(i));%(P)
                Y(s,i,j,t) = q_(i)*norm(f,q(i));%(6)
            end
            Y(s,q0,j,t) = exp(q0_*sum(log(f)));%(6)sc(8)mc
        end
    end
end
% matlabpool close%(8)mc

function [Y,x,int,q,q0,q_1,qix,channels,csel,tinte,tsegs,with1,lx1,ix1,ix2] = mfa_args(x,q,int,til,csel)%(all;)
x = shiftdim(x);%(6.2)(0;4)(1;4.2)(2;4.2)(3:6;)
[lx1,channels] = size(x);
if int(1)>int(end) || int(end)>til || til>lx1 error('In input arguments'), end%(6;)
tinte = til;%(all;)
int = unique(int);%(7.2)
q0 = q==0;
qix = find(q~=0);%(0;4)(1;5)(2;6)(3;7)(4;5)
q = shiftdim(q);%(1;5)(2;6)(3;7)(4;5)
q_1 = 1./q;%(3;7.2)(4;5.2)(5;5.2)(6;)
ix1 = 1:int(end);%(4;1:5,7)
ix2 = [int(end):tinte,ones(1,int(end)-int(1))];%(4;7)
tsegs = floor(lx1/tinte);%(0;4)(1;5)(2;5)(3;6)(4;D,Dt,P,Pt)
%tsegs = lx1-til+1;%(6;)
with1 = tsegs*channels;%(0;4)(1;5)(2,3;)
[i,j] = find(csel);%(C;C)
if true(1)
    csel = [i,j]; csel(:,:,2) = ones(size(csel));%(C;C4)
else
    csel = ones(length(i),2); csel(:,:,2) = [i,j];%(C;C4)
end
Y = zeros(length(int),length(q),channels,tsegs);%(0;4)(3;7)(4;5:7)(5;) @mdfa %
x = reshape(x(1:tinte*tsegs,:),tinte,tsegs*channels);%(4;D,Dt,P,Pt)(5;)
x = cumsum(x-repmat(mean(x),tinte,1));%(0;4)(1;5)(2;5)(3;6)(4;D,Dt,P,Pt)(5;)
