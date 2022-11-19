function Partfn = wtmm_lm = (wtc,q)

 % wtmm_ml implements Wavelet Transform Modilus Maxima method with localmax subfunction
 % wtc are coefitients from wavelet transform.
 % q is the index variable can take any real variable, fot example, q=[-10:1:10].
 % Partfn is partition function.
 % % References:
 % Muzy, J., Bacry, E., Arneodo, A.: Multifractal formalism revisited with wavelets. International Journal of Bifurcation and Chaos in Applied Sciences and Engineering 04 (April 1994)

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

function [Partfn,a1] = wtmm2(wtc,q)
wtc = abs(wtc);a1 = 8;%[wtc,a1] = wtmm_thr(wtc,0.4,8);
maxmap = localmax(wtc,1,false);%true
Partfn = zeros(size(wtc,1),length(q));
for scale_a = size(wtc,1):-1:2
    if ~isempty(find(maxmap(scale_a,:),1))%if nnz(maxmap(scale_a,:)) > 0
        maxln = localmax_(maxmap,scale_a);
        maxwtc = zeros(scale_a,nnz(maxln(1,:)));%not obligatory
        ptr_next = find(maxln(1,:));
        supremum = zeros(1,length(ptr_next));
        for i = 1:scale_a
            ptr_next = maxln(i,ptr_next);%ptr_next(i+1,:)=maxln(i,ptr_next(i,:))
            ptr_corr = ptr_next - 1;%correct position: shl(maxln) - 1
            maxwtc(i,:) = wtc(i,ptr_corr);%maxline(s)WaveletTransformCoefitients
            supremum = max(maxwtc(i,:),supremum);%supremum correction
            maxwtc(i,:) = supremum;
        end
        maxmap(maxln~=0) = 0;%(2)%I=find(maxln);maxmap(I)=0;%erase maxima line(s)(1)
        Maxwtc = repmat(maxwtc,[1 1,length(q)]);
        Q = repmat(shiftdim(q,-1),size(maxwtc));
        partfn = squeeze(sum(Maxwtc.^Q,2));%partfn(a,q)
        Partfn(1:size(partfn,1),:) = Partfn(1:size(partfn,1),:) + partfn;
    end
end

%LOCALMAX Compute local maxima.
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 05-Oct-96.
%   Last Revision: 20-Dec-1999.
%   Copyright 1995-2002 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2002/03/28 17:26:22 $
function x = localmax_(x,rInit)
%function x = localmaxchn(x,rInit)
% Chain maxima - Eliminate "false" maxima.
%-----------------------------------------
ideb = rInit ; step = -1; ifin = 1;
max_down = find(x(ideb,:));
x(ideb,max_down) = max_down;
if rInit<2 , return; end
for i = ideb+step:step:ifin
    max_curr = find(x(i,:));
    val_max  = zeros(size(max_curr));
    for k = 1:length(max_down)
        [nul,ind] = min(abs(max_curr-max_down(k)));
        val_max(ind) = max_down(k);
    end
    x(i,max_curr) = val_max;
    %max_down = max_curr(find(val_max));%(1)
    max_down = max_curr(val_max~=0);%(2)
end
