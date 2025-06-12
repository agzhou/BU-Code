clear all;
files0 = dir('F:\0320_2021_BL2\CP2\*.mat');
files = tsort(files0);
len = length(files);

for i = 1:len
    oldname = files(i).name;
    Oldname = [files(i).folder,'\',files(i).name];
    STR = strsplit(oldname,'-');
    Oldname0 = ['"',Oldname,'"'];
    newname = strcat(strjoin({STR{1:7},num2str(i)},'-'),'.mat');
    command = ['rename' 32 Oldname0 32 newname];
    status = dos(command);
    if status == 0
        disp([oldname, ' rename ', newname]);
    else
        disp([oldname, ' rename failed! ']);
    end
end

function files = tsort(files0)
[~, ind]=sort([files0(:).datenum],'ascend');
files = files0(ind);
end