function mfa(method,s_dir,d_dir,q,int,til,prec,filt,csel,dim)%(3)

 % mfa - Multifractal analysis of time series
 %
 % The parameter method defines multifractal analysis (detrending) method: 'mdfa','mdfast','mdcca','wtmm' or 'prmes'. 'prmes' calculates multifractal parameters.
 % s_dir is the source directory where the structures of data (time series) are. d_dir is the destination directory where the results from the analysis are stored.
 % q is the index variable, for example, q=[-10:1:10].
 % int is the detrending interval, for example, int=round(logspace(log10(10),log10(length(x)/4),32)).
 % til is the width of the time sliding window.
 % prec the is precision, single or double.
 % The parameter filt defines possible filtering of the time series, if it is 'false', then there is no filtering, otherwise filt is a structure with the following fields: filt.n - filter name: 'besel','butter',... ; filt.o - filter order; filt.f - cutoff frequency; filt.s - sampling frequency;
 % When 'mdfast' is used, csel selects which rows of the matrix will be crosscorelated. For example, csel=[1 0; 0 1] means that only the first and the fourth rows are crosscorelated.
 % The parameter dim defines an additional preprocessing of the time series; if dim(1) is true or 1, then an analytical signal is derived from the time series; if dim(2) is true or 1, then a random permutation (shuffling) is performed.
 %
 % Example:
 % Folders and subfolders can be obtained from https://www.physionet.org/content/nesfdb/1.0.0/ -> https://www.physionet.org/content/nesfdb/1.0.0/nesfdb.tar.gz
 % They are placed in folder /work/src/. Then three steps are processed:
 % 1. mfa('mdca','work/src','work/dst',[-16:16],round(logspace(log10(10),log10(450),32)),1800,'single',false,[0 1; 0 0],false(2))
 % 2. fileslink([dst,'/pose/'],'work/om/',{'eh','yng'},{'NULL','STIM'})
 % 3. mfa('prmes','work/om/mdca','work/mfa/pose/mdca',[-16:16],detrint(10,450,32),1800,'single','mdca',[0 1; 0 0],detrint(10,450,32))

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

file=dir(s_dir); file=file(3:end);
%Dat=load(fullfile(s_dir,file(1).name)); Var=fieldnames(Dat); sizx=size(Dat.(Var{1}));%dim=3;%(1)
% X=load(fullfile(s_dir,file(1).name)); sizx=size(X.(cell2mat(fieldnames(X))));%(2)
% fitint=zeros(length(file),sizx(dim),2,'uint16'); h2=zeros(length(file),sizx(dim),'single'); sw=h2;%(1) =
%Fitint=zeros(length(file),sizx(dim),2,'uint16');
%H2=zeros(length(file),sizx(dim),'single'); Sw=H2;%(2)
for i = 1:length(file)
    fs=fullfile(s_dir,file(i).name); [p,fn]=fileparts(fs); fd=fullfile(d_dir,fn); disp(fd)%(5)
    if file(i).isdir%(5)
        if ~strcmp(d_dir,s_dir), mkdir(d_dir,file(i).name), end%(5)
        mfa(method,fs,fd,q,int,til,prec,filt,csel,dim)%(5)
    else%(5)
        %%Dat=load(fullfile(s_dir,file(i).name)); Var=fieldnames(Dat); X=Dat.(Var{1})%(1)
        %% X=load(fullfile(s_dir,file(i).name)); X=X.(cell2mat(fieldnames(X)));%(2)
        %X=load(fs); if isstruct(X), X=X.(cell2mat(fieldnames(X))); end%(5)
        %%if ndims(X)>2, X=X(:,:,logical(csel(1:size(X,ndims(X)-1),:))); end%(4) = err
        %[zi,zj]=find(csel,1,'last'); if ndims(X)>2, X=X(:,:,logical(csel(1:zi,1:zj))); end%(5.2) =
        X = mfa_pref(fs,csel);%(7.2) [X,int]=mfa_pref(fs,csel);%(9)
        switch method%(3)
            case {'mdfa','mdfast','mdcca','wtmm'}%(9) {'mdfa','mdcca','wtmm'}%(7) {'mdfa','mdcca'}
               %[Y,til_] = mfa_prep(X,prec,filt);%() -
                [Y,til_] = mfa_prep(X,prec,filt,dim);%(6) -
                if isinf(evalin('base','til')),til=til_;int=round(logspace(log10(10),log10(til/4),32));end%(9)strcmp(til,'a'),til=='a' -
                [F,    ] = mfa_dtrend(method,Y,q,int,til,prec,csel);%(8)[F,intf]=%() -
                %     fitint(i,:,:) = mfa_fitsel(X,int,dim,file(i).name,q);%(1,7) =
                %%fitint=mfa_fitsel(X,int,dim); Fitint(i,:,:)=fitint;%(2)
                %%if isstruct(filt), fitint(i,:,:)=mfa_fitsel(X,int,dim,file(i).name,q);%(5.2) =
                %%else fitint(i,:,:)=repmat(filt([1,end]),size(X,dim),1); end%(5.2) =
                %      [h2(i,:),sw(i,:)] = mfa_prms(filt,X,q,int,fitint(i,:,:),dim);%(7)(X,q,int,fitint(i,:,:),dim);%(1) =
                %%[h2,sw]=mfa_prms(X,q,int,fitint,dim); H2(i,:)=h2; Sw(i,:)=sw;%(2)
                % [Y(i,:),Z(i,:)] = mfa_statan(h2,sw,int,fitint);%()
                % [Y] = mfa_post(X,varargin);%()
                save(fd,'F')%(5) save(fullfile(d_dir,file(i).name),'F')%(1,2) -
                %%     save(fd,'fitint','h2','sw')%(5) save(fullfile(d_dir,file(i).name),'fitint','h2','sw')%(1,2) =
            case 'prmes'
                %%[Y,til ] = mfa_prep(X,prec,filt);%() -
                %     [Y,til_] = mfa_prep(X,prec,filt,dim);%(6) -
                %     if isinf(evalin('base','til')),til=til_;int=round(logspace(log10(10),log10(til/4),32));end%(9)strcmp(til,'a'),til=='a' -
                %     [F,    ] = mfa_dtrend(method,Y,q,int,til,prec,csel);%(8)[F,intf]=%() -
                fitint(i,:,:) = mfa_fitsel(X,int,dim,file(i).name,q);%(1,7) =
                %fitint=mfa_fitsel(X,int,dim); Fitint(i,:,:)=fitint;%(2)
                %if isstruct(filt), fitint(i,:,:)=mfa_fitsel(X,int,dim,file(i).name,q);%(5.2) =
                %else fitint(i,:,:)=repmat(filt([1,end]),size(X,dim),1); end%(5.2) =
                [h2(i,:),sw(i,:)] = mfa_prms(filt,X,q,int,fitint(i,:,:),dim);%(7)(X,q,int,fitint(i,:,:),dim);%(1) =
                %[h2,sw]=mfa_prms(X,q,int,fitint,dim); H2(i,:)=h2; Sw(i,:)=sw;%(2)
                % [Y(i,:),Z(i,:)] = mfa_statan(h2,sw,int,fitint);%()
                % [Y] = mfa_post(X,varargin);%()
                %     save(fd,'F')%(5) save(fullfile(d_dir,file(i).name),'F')%(1,2) -
               % save(fd,'fitint','h2','sw')%(5) save(fullfile(d_dir,file(i).name),'fitint','h2','sw')%(1,2) =
        end
    end%(5)
end
% save(fullfile(d_dir,'h2sw'),'fitint','h2','sw')%(1) save(fullfile(d_dir,'h2sw'),'Fitint','H2','Sw')%(2) =
if exist('h2','var'), save([d_dir,'/',fn([1:end-9,end-4:end])],'fitint','h2','sw'), end%(5.2) save(fd,'fitint','h2','sw')%(5) =

function X = mfa_pref(fs,csel)%(7.2),'201306061209.00'

% Subfunction to mfa. Conversion of the input data.

%Dat=load(fullfile(s_dir,file(i).name)); Var=fieldnames(Dat); X=Dat.(Var{1})%(1)
% X=load(fullfile(s_dir,file(i).name)); X=X.(cell2mat(fieldnames(X)));%(2)
X=load(fs); if isstruct(X), X=X.(cell2mat(fieldnames(X))); end%(5)
%if ndims(X)>2, X=X(:,:,logical(csel(1:size(X,ndims(X)-1),:))); end%(4) = err
[zi,zj]=find(csel,1,'last'); if ndims(X)>2, X=X(:,:,logical(csel(1:zi,1:zj))); end%(5.2) =
%X=X(16384:287228,1:end-1);%eyes+head

function [Y,lengthy] = mfa_prep(X,prec,filt,dim)%(4)

% Subfunction to mfa. Pre-processing and transformation: filters, extracting the amplitude and phase, surrogate data.

% function [Y,lengthy]=mfa_prep(X,prec,filt)%(3)
if nargin<4, dim=false(2); end, if nargin<3, filt=false; end, if nargin<2, prec='double'; end%(4.)
precx = class(X); X = cast(X,prec);%(3) X=double(X);
%[b,a]=butter(5,10/(100/2));%(2) %X=filt_m(X,'butter',5,10,100,'low');%(1)
% [b,a]=feval(filt.n,filt.o,filt.f/(filt.s/2)); X=filtfilt(b,a,X);%(3)%(2:3)
if isstruct(filt), [b,a]=feval(filt.n,filt.o,filt.f/(filt.s/2)); X=filtfilt(b,a,X); end%(4)201306020816.34
if exist('dim','var')&&dim(1), X=[X,analysig(X)]; end%(4) X=[X,analysig(X)];
if exist('dim','var')&&length(dim)==2&&dim(2), X=[X,rnd_all('perm',X,1)]; end%(4) X=[X,rnd_all('perm',X,1)];
Y = cast(X,precx);%(3) Y=single(X);
lengthy = size(Y,1);

function varargout = mfa_dtrend(method,D,q,int,til,prec,csel)%(4)201306100944.1

% Subfunction to mfa. Uniform interface to core functions for multifractal analysis that achieves seamless addition of new methods.

% function [F,intf] = mfa_dtrend(method,D,q,int,til,prec,csel)%(,varargin)%()'-fwv
%[F,intf] = mdfa(D,q,int,til,prec);%(mdfa)%(1)mdfa8([D,analysig(D)],q,int,length(D),ones(2));%(9.6.fast)
%[F,intf] = mdcca(D,q,int,til,prec,csel);%(mdcca)%(1)
% [F,intf] = feval(method,D,q,int,til,prec,csel);%(,varargin{1})%(2)'-fwv
switch nargout%(4){"=
    case 1
        varargout(1) = {feval(method,D,q,int,til,prec,csel)};
    case 2
        [F,int]=feval(method,D,q,int,til,prec,csel); varargout(1)={F}; varargout(2)={int};%(2,4)(4)(4)
end%(4)}"=
% F = F(:,:,logical(csel));%(,varargin{1})%(3)
%F=shiftdim(F,4); F=F(:,:,:,logical(csel)); F=shiftdim(F,1);%(csel,t)
%F=shiftdim(F,ndims(F)-5); F=F(:,:,:,logical(csel)); F=shiftdim(F);%(t,csel)

%'~/My Documents/My Work/work/dfa/m/fit_sel.m','14.04.08','./fitsel.m','19.03.07,12:28'
%'~/My Documents/My Work/work/dfa/m/mfa_fitsel.m',
%'(2)01.11.10,22:31; (3)02.11.10,03:01,04:12; 03.11.10,02:17; (4)04.04.11,20:49; (5)08.06.11,21:20;...
% (5.2)09.06.11,21:44; (6)21.02.12,06:45; (7)201306061110.05 (8)201306070409.0'
function Fitint = mfa_fitsel(X,int,dim,filename,q)%(5.2)

% Subfunction to mfa. Linear approximation of Fq(s) (Zq(s)) - automatically or with the GUI selection of the range for obtaining a faultless multifractal spectrum.

%function Fitint=mfa_fitsel(s_dir,d_dir,int)%(2) =fitsel(s_dir,d_dir,int)%(1)
%file = dir(s_dir);file = file(3:end);%(1,2)
%Dat=load(fullfile(s_dir,file(1).name));Var=fieldnames(Dat);%(2) load(fullfile(s_dir,file(1).name));%(1)
%[zi,zj]=find(csel,1,'last'); if ndims(X)>2, X=X(:,:,logical(csel(1:zi,1:zj))); end%(5.2)% X = X(:,:,logical(csel));%(4)
if isscalar(dim)%(7)
Fitint=zeros(size(X,dim),2,'single');%(6) =zeros(size(X,dim),2);%(3) =zeros(length(file),size(Dat.(Var{1}),3),2);%(2)...
    % =zeros(2,size(D,2),length(file));%(1)
%for i = 1:length(file)%(1,2)
    %Dat=load(fullfile(s_dir,file(1).name));Var=fieldnames(Dat);%(2) load(fullfile(s_dir,file(i).name));%(1)
    for j = 1:size(X,dim)%(3) =1:size(Dat.(Var{1}),3)%(2) =1:size(D,2)%(1),disp(file(i).name),disp(j)%(1,2)
		button = 'No';%(5)
		while strcmp(button,'No')%(5)
        loglog(int,X(:,:,j)),title([filename,num2str(j)])%(3) loglog(int,Dat.(Var{1})(:,:,j))%(2) loglog(int,squeeze(D(:,j,:)))%(1)
        [fitint,y] = ginput(1); fitint(fitint<int(1)) = int(1);%(2) [Fitint(:,j,i),y]=ginput(2)%(1)
        loglog(int(int>=fitint),X(int>=fitint,:,j))%(3) loglog(int(int>=fitint),Dat.(Var{1})(int>=fitint:end,:,j))%(2)
        [fitint(2),y] = ginput(1); fitint(fitint>int(end)) = int(end); fitint(fitint<int(1)) = int(1); %(2) %(6)
		fitint = sort(floor(fitint)); intf = int>=fitint(1) & int<=fitint(end);...
            [a,f] = spect(tau(fit(X(intf,:,j),int(intf),'log',1),q),q); plot(a,f)%(5.2)
		button = questdlg('Next?'); assert(~strcmp(button,'Cancel'),'Cancel')%(5)
		end%(5)
		fitint(fitint==fliplr(fitint)) = NaN; Fitint(j,:) = fitint;%(6) =single(sort(floor(fitint)));%(3)...
            % F..(i,j,:)=single(sort(floor(fitint)));%(2))
    end%,%fitint=single(sort(floor(Fitint(:,:,i))));save(fullfile(d_dir,file(i).name),'fitint')%(1)
    %index = int>=fitint(1) & int<=fitint(end);%(3)
%end%(1,2),%Fitint = sort(floor(Fitint));close%close(gcf)%(1)
else Fitint=dim([1,end]);%(8)=repmat(dim([1,end]),size(X,3),1);%(7)
end%(7)


function [h2,sw] = mfa_prms(method,X,q,int,fitint,dim)%(5)

% Subfunction to mfa. Function for obtaining multifractal parameters.

if ~isscalar(dim), dim=3; end%(7)
h2=zeros(1,size(X,dim),'single'); sw=h2; fitint=shiftdim(fitint,1);%(4)fitint)%(3) zeros(1,size(X,dim))%(2)zeros(1,size(fitint,1))%(1)
% for i = 1:size(X,dim)% size(fitint,1)%()
%     switch method%(5)
%         case {'mdfa','mdcca'}%(5)
%             intf = int>=fitint(i,1) & int<=fitint(i,end);%()
%             hq = fit(X(intf,:,i),int(intf),'log',1); hq(hq==0) = NaN; %() %(3)
%             tauq = tau(hq,q);%()
%         case {'wtmm'}%(5)
%             tauq=fit(X(:,:,i),int,'lin',1); hq=(tauq+1)./repmat(q(:),1,size(tauq,2));%(5)
%     end%(5)
%     h2(i) = hq(q==2);%() h2(i,:)=hq(q==2,:);%()
%     a = spect(tauq,q); a = a(isfinite(a)); %()' %(3.2)"
%     sw(i)=a(1)-a(end);%() sw(i,:)=a(1,:)-a(end,:);sw1=sw(i,:);%() =a(find(isfinite(a),1))-a(end);%(3.2)
% end
switch method%(6)
    case {'mdfa','mdfast','mdcca'}%(11) {'mdfa','mdcca'}%(5)
        if ~isvector(fitint)%(8)
            hq = zeros(size(X,dim-1),size(X,dim),'single');%(6)
            for i = 1:size(X,dim)% size(fitint,1)%()
                intf = int>=fitint(i,1) & int<=fitint(i,end);%()
                hq(:,i) = fit(X(intf,:,i),int(intf),'log',1);% hq(hq==0) = NaN; %() %(3)
            end
        else intf=int>=fitint(1)&int<=fitint(end); hq=fit(X(intf,:,:,:),int(intf),'log',1);%(8)
        end%(8)
        hq(hq==0) = NaN; tauq = tau(hq,q); %(8)' %()"
    case {'wtmm'}%(5)
        if size(X,2)~=length(q)||size(X,1)~=length(int), int=X(:,1,1); X=X(:,2:end,:); q=q(1):0.2:q(end); end%(10)
        tauq = fit(X,int,'lin',1); hq = (tauq+1)./repmat(q(:),1,size(tauq,2));%(6)
end%(6)
h2 = hq(q==2,:,:);%(6) h2(i)=hq(q==2);%() h2(i,:)=hq(q==2,:);%()
a = spect(tauq,q); %a = a(isfinite(a)); %()' %(3.2)"
sw = a(1,:,:)-a(end,:,:);%(6,9) sw(i)=a(1)-a(end);%()sw(i,:)=a(1,:)-a(end,:);sw1=sw(i,:);%()=a(find(isfinite(a),1))-a(end);%(3.2)

% function [Y,Z] = mfa_statan(h2,sw,int,fitint,dim);%()
% [Y(i,:),Z(i,:)] = mfa_statan(h2,sw,int,fitint,dim);%()

% function [Y] = mfa_post(X,varargin);%()
% h2(int,:) = h21(int,:); sw(int,:) = sw1(int,:); fitint(int,:,:) = fitint1(int,:,:);
