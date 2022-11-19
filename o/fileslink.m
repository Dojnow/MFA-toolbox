function fileslink(s_dir,d_dir,P,C) %conv_rec5()

 % fileslink tracking down recursively and searching in the tree  s_dir of the file system, selecting and creating links to files in folder d_dir grouped by various indications, such as experimental conditions and test subjects. s_dir is the source directory, d_dir is the destination directory.
 % P and C are the indications.
 % See also an example in the description of the function mfa where P={'eh','yng'}; C={'NULL','STIM'};.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

% file select and link %dst='work.mat/dst';
for p=P%{'eh','yng'}%(cr.2) {'group_1','group_2','group_3'}%(cr.4)
    for c=C%{'NULL','STIM'}%(cr.2) {'c','f','m','s'}%(cr.4)
        %conv_rec4([dst,'/pose/',p{1}],['work.mat/om/',p{1},',',c{1}],c{1})%(0)
        conv_rec4([s_dir,p{1}],[d_dir,p{1},',',c{1}],c{1})%(1)
    end
end

function conv_rec4(s_dir,d_dir,c_dir)%recsec
% c_dir - files in dir.s with name c_dir, which are within s_dir, are copied to d_dir
if ~exist(d_dir,'dir'), mkdir(d_dir), end%(3)
file=dir(s_dir); file=file(3:end);
for i = 1:length(file)
    fs=fullfile(s_dir,file(i).name);% disp(fs)
    [p,casedir]=fileparts(s_dir);%(4)
    if file(i).isdir% && ~strcmp(file(i).name,c_dir)%(1.1)
        %if strcmp(file(i).name,c_dir), conv_rec4(fs,d_dir,fs), else conv_rec4(fs,d_dir,c_dir), end%(1.2)
        conv_rec4(fs,d_dir,c_dir)%(1.1,2)
    %elseif file(i).isdir && strcmp(file(i).name,c_dir), conv_rec4(fs,d_dir,fs)%(1.1)
    %elseif isempty(strfind(fs,c_dir)), return%(2)
    elseif ~strcmp(casedir,c_dir), return%(4)
    else%(2)
        unix(['ln -s ~/My\ Documents/My\ Work/work/"',fs,'" "',d_dir,'"']);%(3) %copyfile(fs,d_dir)%(1,2)
    end%(1,2)
end
