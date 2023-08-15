function beat_excel(xlpath,topwhitespace)
%
%   beat_excel(xlpath,topwhitespace)
%
% Function takes the resuits*.txt produced from the python analysis of
% mechanical loss meaurements and writes them to a specific excel sheet.
% Function also generates excel formulae to calculate averages and standard
% deviations of  mechanical loss values for each frequency. This relies on
% the k-means clustering algorythim is required.
%
% The excel sheet which is created will have the same name as the folder in
% which the results*.txt file is found
%
%
% INPUTS
%   xlpath          String pointing to the excel file which should be
%                   written too.
%
%   topwhitespace   scalar value determining the number of blank cells at
%                   the top of the excel sheet before the data is written.
%                   default = 3
%
%
%
%
%
% Author S.tait 2023
% simon.tait@glasgow.ac.uk
%
% 
% v1.0.1    - Bug fixes and changes to writecell and writematrix functions
%             to enable writing to excel on windows systems. Also added a cleaning up
%             system which zips all images and pickle files after user review to save
%             space.
%

if ~exist('xlpath','var')
    fprintf('\nSpecify which excel file to store outputs\n')
    [tmpfile,xlpath] =uigetfile('*.xlsx','Specify which excel file to store outputs') ;
    xlpath=  fullfile(xlpath,tmpfile) ;
end

% if exist(xlpath) ~= 7
%     cpritnf('err' ,'ERR: Path specified does not point to a folder\n')
%     cpritnf('err' ,'Please make sure xlpath is a string and points to the folder\n')
%     cpritnf('err' ,'which contains the desired excel file')
% end

if iscell(xlpath)
    xlpath = char(path);
end

if ~exist('topwhitespace','var')
    topwhitespace = 3;
end

%

% xlpath =  '/Users/igradmin/Dropbox/Mechanical Loss-BERNERAY/3_inch_Disk_Samples/G&H_Disks_Batch_30070698/G&H_Disk_30070798-04';
% xlfile ='Disk04_TiSiO2HR_6040.xlsx';


%split path into parts
[xlpath,xlfile,extension] = fileparts(xlpath);
xlfile =strcat(xlfile,extension);


dirparts = strsplit(pwd,filesep) ;
sheetname = dirparts{end-1};

% check that sheetname is  not over max name size
if numel(sheetname)>31
    %check if name contains 'suspension'
    if contains(lower(sheetname),'suspension')
        sheetname = strrep(lower(sheetname),'suspension','sus');
    end
    %check length
    if numel(sheetname)>31
        sheetname = string(today('datetime'));
        cprintf('hyper','\n Folder name is too long, sheet will be named %s','sheetname')
    end
end



%find and load results data
results = readtable(char(fnames('results*.txt')));
results = results(:,[1,3:end]);

results = table2array(results);

%load pre-beat / Total analysis losses
try
    pre = load(fullfile(pwd,'PyTotalAnalysis.mat'),'phi');
    
catch
    cprintf('err','\n Could not find ''PyTotalAnalysis.mat''\n')
    [MatFile] =uigetfile('*.mat','Please Specifiy which file to read');
    try
        pre = load(fullfile(pwd,MatFile),'phi');
    catch
        cprintf('err','\nPython outputs could not be loaded ''\n')
        pre  = nan(size(results(:,1)));
    end
    
end
results = [results,pre.phi(results(:,end-2)+1)'];

results =  sortrows(results,1);


%define english alphabet for manovering in execel
alphab =  {'a','b', 'c', 'd', 'e', 'f', 'g','h', 'i', 'j', 'k','l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
start_col = find(cellfun(@(x) any(strcmp({'a'}, x)), alphab));

xlheadings = {'Freq', 'Loss 1','Loss 2','low err 1','up err1','low err 2','up err 2','loop num','good resid','use 2 qs','TotalAnalyis','USER Review'};
xlrange_heading = char(strcat('B',string(topwhitespace-1),':',capitalize(alphab(size(results,2)+2)),string(size(results,1)+topwhitespace)));

% prep data for clustering analysis
freq = results(:,1);

% rough_f =unique(ceil(freq./100)*100);
counts = groupcounts(ceil(freq./100)*100);

%perform clusering analysis and generate groups and labels
labels = kmeans(freq,numel(counts));
groups = unique(labels,'stable');

%sort groups so that frequencies occur in order
% tmp = unique(sortrows([labels,ceil(freq./100)*100],2),'stable');
% groups = tmp(tmp<50);

for j=1:numel(groups)
    counts2(j) = numel(find(ismember(labels,groups(j))));
end
counts = counts2';

%write headings to excel
fprintf('\n Writing Data...')

if ispc ==1 
    writecell(xlheadings,fullfile(xlpath,xlfile),'sheet',sheetname,'Range',xlrange_heading)
    % chk = xlswrite(fullfile(xlpath,xlfile),xlheadings,sheetname,xlrange_heading);
else 
    chk = xlwrite(fullfile(xlpath,xlfile),xlheadings,sheetname,xlrange_heading);
end 



%calculate starts and ends of each frequency for nice excel formatting
ends = cumsum(counts)+2 +[0,1:numel(counts)-1]';
starts = (ends-counts)+1;

starts = [starts(1);starts(2:end)+1];




for i =1:numel(groups)
    groupidx(:,i) = dynamic(find(ismember(labels,groups(i))),max(counts),'NaN');
    excelidx(:,i) = dynamic([starts(i):ends(i)]',max(counts),'NaN');
end


fitcheck = string(nan(max(counts),numel(groups)));
for i =1:numel(groups)
    %find indicies of each frequency in results
    %sort so that negative values are at the bottom
    outdata = results(groupidx(~isnan(groupidx(:,i)),i),:);
    
    %generate names for fit figures produced by python
    % these will be used to hyperlink - ADD THIS LATER USING DROPBOX API
    fitnames = strcat('fit_',string(outdata(:,8)),'.png');
    for k =1:numel(fitnames)
        hypernames{k} = strcat('=HYPERLINK("',which(fitnames{k}),sprintf('","%s")',fitnames{k}) );
        %user checking fit is okay
        try
            imshow(sprintf('%s',which(fitnames{k})),'InitialMagnification',300)
            pause(0.5)
            fitcheck{k,i} = input('\nAre you happy with this fit [Y/N] \n','s');
            
            
            if isempty(fitcheck{k,i}) || contains( lower(fitcheck{k,i}),'y')
                
                fitcheck{k,i} ='1';
            
            else
                
                fitcheck{k,i} ='0';
            end
            
            close
        catch
            fitcheck{k,i} ='0';
            fprintf('\n ERROR file %s could not be loaded\n',char(fitnames{i}))
        end
        
    end
    
    xlrange_names = char(strcat('A',string(starts(i)),':','A',string(ends(i))));
    
    
    
    %write file names to excel
    if ispc
        % chk3 = xlswrite(fullfile(xlpath,xlfile),hypernames',sheetname,xlrange_names);
        writecell(hypernames',fullfile(xlpath,xlfile),'sheet',sheetname,'Range',xlrange_names)

    else
        chk3 = xlwrite(fullfile(xlpath,xlfile),hypernames',sheetname,xlrange_names);
    end


    %generate excel cells
    xlrange_data =char(strcat('B',string(starts(i)),':',capitalize(alphab(size(results,2)+2)),string(ends(i))));
    outdata = [outdata,str2double(fitcheck(1:size(outdata,1),i))>0];
    
    %write data to excel
    if ispc 
    % chk2 = xlswrite(fullfile(xlpath,xlfile),outdata,sheetname,xlrange_data);
    writematrix(outdata,fullfile(xlpath,xlfile),'sheet',sheetname,'Range',xlrange_data)
    else 
    chk2 = xlwrite(fullfile(xlpath,xlfile),outdata,sheetname,xlrange_data);
    end 
    average_idx = excelidx(outdata(:,end)>0,i);
    try
        if ~isempty(average_idx)
            outform{i,1} = strcat('=AVERAGE(',csvstring(strcat('B',string(average_idx))),')');
            pause(0.5)
            outform{i,2} = strcat('=AVERAGE(',csvstring(strcat('C',string(average_idx))),')');
            pause(0.5)
            outform{i,3} = strcat('=AVERAGE(',csvstring(strcat('D',string(average_idx))),')');
            pause(0.5)
            outform{i,4} = strcat('=STDEV(',csvstring(strcat('C',string(average_idx))),')');
            pause(0.5)
            outform{i,5} = strcat('=STDEV(',csvstring(strcat('D',string(average_idx))),')');
            pause(0.5)
            outform{i,6} = strcat('=AVERAGE(',csvstring(strcat('L',string(excelidx(find(~isnan(outdata(:,end-1))),i)))),')');
            pause(0.5)
            outform{i,7} = strcat('=STDEV(',csvstring(strcat('L',string(excelidx(find(~isnan(outdata(:,end-1))),i)))),')');
        else
            outform{i,1} = 'NaN';
            pause(0.5)
            outform{i,2} = 'NaN';
            pause(0.5)
            outform{i,3} = 'NaN';
            pause(0.5)
            outform{i,4} = 'NaN';
            pause(0.5)
            outform{i,5} = 'NaN';
            pause(0.5)
            outform{i,6} = 'NaN';
            pause(0.5)
            outform{i,7} = 'NaN';
        end
    end
    
    
    %construct averages array
    
    clear out2 out3 fitnames outnames outdata average_idx
    
    
    cprintf('.')
    
end

%generate formulae for averaging and STD of mecahnical losses
% for j = 1:numel(starts)
%     out(j,:) = {string(sprintf('=AVERAGE(B%d:B%d)',starts(j),ends(j))),string(sprintf('=AVERAGE(C%d:C%d)',starts(j),ends(j))),string(sprintf('=AVERAGE(D%d:D%d)',starts(j),ends(j))),string(sprintf('=STDEV(C%d:C%d)',starts(j),ends(j))),string(sprintf('=STDEV(D%d:D%d)',starts(j),ends(j)))};
% end

if ispc ==1 
    write_range = char(strcat('N',string(topwhitespace-1),':',capitalize(alphab(find(cellfun(@(x) any(strcmp({'n'}, x)),alphab))+size(outform,2))),string(size(outform,1)+topwhitespace-1)));
    table_headings = {'Freq', 'PyLoss 1','PyLoss 2','STDEV PyLoss1','STDEV PyLoss2','ToA Loss','STDEV ToA'};

    writecell(table_headings,fullfile(xlpath,xlfile),'sheet',sheetname,'Range',write_range)


    % xlswrite(fullfile(xlpath,xlfile),table_headings,sheetname,write_range);
    pause(0.5)
else 
    write_range = char(strcat('N',string(topwhitespace-1),':',capitalize(alphab(find(cellfun(@(x) any(strcmp({'n'}, x)),alphab))+size(outform,2))),string(size(outform,1)+topwhitespace-1)));
    table_headings = {'Freq', 'PyLoss 1','PyLoss 2','STDEV PyLoss1','STDEV PyLoss2','ToA Loss','STDEV ToA'};

    xlwrite(fullfile(xlpath,xlfile),table_headings,sheetname,write_range);


    pause(0.5)
end
%write formulae to excel
if any(sum((cellfun(@isempty,outform)))>1)
    outform(find(cellfun(@isempty,outform))) = {'NaN'};
end

if ispc ==1 
    
    outform_range = char(strcat('N',string(topwhitespace),':',capitalize(alphab(find(cellfun(@(x) any(strcmp({'n'}, x)),alphab))+size(outform,2))),string(size(outform,1)+topwhitespace-1)));
    
    writecell(outform,fullfile(xlpath,xlfile),'sheet',sheetname,'Range',outform_range)
    % xlswrite(fullfile(xlpath,xlfile),outform,sheetname,outform_range);


    pause(0.5)
else 
    outform_range = char(strcat('N',string(topwhitespace),':',capitalize(alphab(find(cellfun(@(x) any(strcmp({'n'}, x)),alphab))+size(outform,2))),string(size(outform,1)+topwhitespace-1)));
    xlwrite(fullfile(xlpath,xlfile),outform,sheetname,outform_range);


    pause(0.5)
end 



cprintf('text','Done\n')

pause(5) 

reply =input('Do you wish to open the excel file ? [Y/N]','s');
if contains(lower(reply),'y')
system(sprintf('! open %s', char(fullfile(xlpath,xlfile))    ))
end 


fprintf('\n Cleaning up directory...\n')
% Step 1: Find the files
filePattern = fullfile('.', 'fit*.png');
fileList = dir(filePattern);

% Step 2: Create the zip archive
zipFilename = 'fit_images.zip';
filesToZip = strcat({fileList.folder}, filesep, {fileList.name});
zip(zipFilename, filesToZip);

% Step 3: Delete the images
for i = 1:numel(fileList)
    delete(fullfile(fileList(i).folder, fileList(i).name));
end

% Step 1: Find the files
filePattern = fullfile('.', '*.pickle');
fileList = dir(filePattern);

% Step 2: Create the zip archive
zipFilename = 'Pickles.zip';
filesToZip = strcat({fileList.folder}, filesep, {fileList.name});
zip(zipFilename, filesToZip);

% Step 3: Delete the images
for i = 1:numel(fileList)
    delete(fullfile(fileList(i).folder, fileList(i).name));
end







end
