function  conv_dir(s_dir,d_dir)

 % Converts an ASCII.dat files to a mat files format.
 %
 % s_dir is the source directory, d_dir is the destination directory.

 % Copyright (C) 2020 Peter Dojnow
 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

file = dir(s_dir);%(1:7)file = file(3:end);%(8,10:11,15,16.2)
for i = 3:length(file)%(3:4,7,9)=1:length(file)%(1:2,8,10:11,15,16.2) =1:9%(14)
    fs = fullfile(s_dir,file(i).name);%(1:7,9:11,15,16.2)
    [p,fn,ext] = fileparts(fs);%(1:6,8.2:9,11,13) fn=['noise0',num2str(i)];%(14)
    fd = fullfile(d_dir,fn);%(1:6,8.2:9,11,13)
    D = load(fs); D = single(D(1:end-5,:));%(11.7)
    dat(i).data=D; dat(i).date=file(i).date; [p,dat(i).name]=fileparts(file(i).name);...
    save(fd,'D');%(7:8,11.1:5,13)
    unix(['touch -mt ',datestr(file(i).datenum,'yyyymmddHHMM.SS'),' ',[fd,'.mat']])%(11.6.3)
end
