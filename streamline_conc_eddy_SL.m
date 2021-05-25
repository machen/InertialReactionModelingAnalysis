close all;
clear all; 
clc; 


caseDir = [pwd '\COMSOL data\'];

%% read data
 fid = fopen([caseDir '3d_streamline_re_300.txt'],'r'); %%open file
    
 tmp = fgets(fid); %% read 1st line
while ~contains(tmp,'Nodes') %%moving the cursor to line in which "Nodes" appears 
    tmp = fgets(fid); %%save line to tmp
end
note = textscan(tmp, '%s','delimiter',':'); %%retrieve the line with delimiter ":" to note
nLine = str2num(note{1}{2}); %%retrieve # of streamlines to nLine

while ~contains(tmp, '% x') %%read until "% x" appears 
    tmp = fgets(fid); %%take line containing "% x"
end


data = zeros(nLine,5); %%make a dataset with zeros of nLine rows and 5 columns

for iLine= 1:nLine  %% loop from iLine from 1 to nLine(total streamline number)  
    s = fgets(fid); %%put cursor to beginning of next line
    data(iLine,:) = sscanf(s, '%f %f %f %d %f',[1 5]); %%take coordinates and data
end

fclose(fid);
%% plot data
clf
eddyLineNo = []; %%make empty variable
deanflowLineNo = []; %% make empty variable
indLine = unique(data(:,4)); %%retrieve unique number for streamlines (from column #4 of 'data' and list that into single column
for iLine = 1:numel(indLine) %%loop from 1 to # of streamlines
    x = data((data(:,4)==indLine(iLine)),1); 
    y = data((data(:,4)==indLine(iLine)),2); 
    z = data((data(:,4)==indLine(iLine)),3); 
    c = data((data(:,4)==indLine(iLine)),5); 
    if sum(c>0.01)>0 
        if sum(y>50)>0 %%if y contains at least one value that is greater than 50
            yPart = y(y>50); %%save values greater than 50 to yPart
            tmp = diff(yPart); %%
            if sum(tmp<0)>0
                eddyLineNo = [eddyLineNo; iLine];
                plot3(x,y,z,'r','lineWidth',1.5); hold on
                xlim([-1 1]*1000)
                ylim([-1 1]*1000)
                zlim([0 1]*70)
                drawnow;
            else
                deanflowLineNo = [deanflowLineNo; iLine];
                plot3(x,y,z,'b','lineWidth',1.5); hold on
                xlim([-1 1]*1000)
                ylim([-1 1]*1000)
                zlim([0 1]*70)
                drawnow;
                
            end
        end
    end
end
