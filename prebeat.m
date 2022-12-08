function prebeat(talk,samplename,auto,overtemp)
%
% Current Version 2.0
%
% This function takes ringdown data from mechanical loss systems and allows
% for logical inspection of the data. A simple 'pass' / 'reject' can be
% invoked allowing data to be kept for further analysis or flat out
% rejected. If there are artifcacts which would cause influence an
% incorrect fit such as long portions of data where the ringdown has
% finished, these can be cropped out. Only the cropped portions of these
% regions are saved for further analysis
%
% On the event of a crash or error the outputs of the program will be saved
% to the working directory. Upon re-running, the function will look
% in the current directory and subfolders for the saved files and
% automatically resume at the crashed point.
%
%
%           prebeat(verbose,'samplename')
%
%  Example:
%           prebeat(1,'TestSample') ;
%
% INPUTS
%
%  verbose         : Turn on prints to command line 1 = on , 0 =off
%
%  samplename      : Name of Sample / Set of data - used for save files and
%                     plots. should be entered as a string
% OUTPUTS
%
% PyTotalAnalysis.mat           :matlab results file which contains the
%                                 final analysed data
%
% Processed_Mode_X.txt          :file  written to a .txt file in the
%                                selected directory containing a breakdown
%                                of measurements for each mechanical  mode
%                                 format = 'freq' 'offset' 'amplitude' 'fitted tau' ' Mechanical loss' 'Origin filename'
%
% Log_<samplenane>_.txt         :Log file containing the results of each
%                                ringdown after it has been presented to the user.
%
% Author  S. Tait 2022
%
% Email   simon.tait@glasgow.gla.ac.uk
%
% v2.0      - Changes to Log file saving structure. Now saves in easy to read
%             format with readtable.m
%           - All outputs from program will be saved on crash or error which
%             can be resumed on re-running prebeat.
% v1.9.3    - Removal of old/new input in favour of automatic directory search
% v1.9.2    - Improved handling of windowed data for python analysis
% v1.9.1    - Improved robustness to code running on different operating
%             systems
% v1.9      - Added filtering to ringdown/scan data to remove large spikes
% v1.8s     - Option to add automatic Ringdown analysis -ALPHA
% v1.8      - Now extracts temperature information from log files
% v1.7      - Removal of Automated Rejection ( Requires tweaking for glitches)
% v1.6      - Recording of fitted loss valvues and record of fitted coordinated
% v1.5      - Auto-rejection of 'just-noise' ringdowns
% v1.4      - Storage of variables optimised for RAM management
% v1.3      - Functionallity on windows systems has been tested and improved.
% v1.2      - Fixed known bug for final output summary on command line
%           - Tweaking of script for diffeent operating systems
%           - ModeNumber functionallity removed
%
currentVersion = 2.0;

% for a = 1:nargin-3  % only process varargin, no other parse_args() arguments
%     thisArg = varargin{a};
%     if ischar(thisArg) && ~isempty(thisArg)
%         switch
%             case '-update'
%                 updateInstalledVersion();
%         end
%     end
% end

%% Check for new version - update if newer version exists
if checkForNewerVersion_MATLAB(currentVersion)
    updateInstalledVersion()
end

%% Look for failed versions in current directory - prompt user to load if exist
loadflag = 0;
% look for recovery files in current directory
currentdirs = dir;
% find all subfolders
subdirs = {currentdirs(find([currentdirs.isdir])).name};
subdirs = subdirs(~contains(subdirs,'.'));
fileidx = isfile(fullfile(subdirs,'Recovered_outputs.mat'));

if sum(isfile(fullfile(subdirs,'Recovered_outputs.mat')))  ==1
    fprintf('\n Failed outputs from a previous run have been found in\n%s\n ',char(strcat('... \.',subdirs(fileidx))))
    reply = input('Do you wish to load this file? [Y\N]','s');
    if contains(lower(reply),'y') || isempty(reply)
        fprintf('Loading previous outputs...\n')
        load(char(fullfile(pwd, subdirs(fileidx),'Recovered_outputs.mat')))
        loadflag=1;

    else
        clc
    end
end

clc

%% Set up log file
if loadflag ==0
    if strcmpi(samplename,'test')
        mkdir Test
        FolderName = fullfile(pwd,'Test');
    else
        FolderName = uigetdir;
    end

    %get operator ID for log files
    [ret, name] = system('hostname');
    if ret ~= 0
        if ispc
            name = getenv('COMPUTERNAME');
        else
            name = getenv('HOSTNAME');
        end
    end
    ComputerID = strtrim(lower(name));
    clear name ret


    %% search current directory for ringdown/scan files and match appropriately
    if numel(fnames('scan*'))>1
        prefix ={'scan_'};
        oldnew =2;
    elseif numel(fnames('spec'))>1
        prefix = {'spec_'};
        oldnew =1;
    else
        error('No files of the ''spec'' or ''scan'' formats could be found Please make sure this script has been run in the correct directory')
        return
    end

    spec_files = dir([char(prefix),'*.txt']);


    ring_files = dir('ring_*.txt');
    ring_files1 = dir('ring_*.txt');

    for i = 1:length(ring_files)
        ring_files(i).name = strrep(ring_files(i).name,regexp(ring_files(i).name,'ring[_]','Match'),'');
        ring_name(i,:) = cellstr(ring_files(i).name);
    end
    clear  i

    for i = 1:length(spec_files)
        spec_files(i).name = strrep(spec_files(i).name,regexp(spec_files(i).name,char(prefix),'Match'),'');
        spec_name(i,:) = cellstr(spec_files(i).name);
    end
    clear i

    %remove files which do not havea successful ringdown
    measured = ismember(spec_name,ring_name);
    spec_name = spec_name(measured == 1 );
    clear tmp

    %initialise output matricies
    ring_times  = [];
    ring_amps  = [];
    ia =[];
    temps = [];

    %find size of largest file  if taking too long to load data - set max
    %length as 100000 - this speeds up analysis when pulling data from
    %cloud services like Dropbox or OneDrive
    tic
    for ii = 1:length(spec_name)

        if toc>60 && length(spec_name)<400
            numLines(ii) = getLines(strcat(char(strcat('ring_',ring_name{ii}))));
        else
            numLines = ones(size(spec_name)).*100000;
            ii =length(spec_name);
        end
    end

    C=0;
    matchk = max(numLines) ; clear ii;
    counter = 0;

    histfile =  fullfile(FolderName,strcat('Log_',samplename,'_',datestr(now,'DD-mm-YYYY'),'.txt'));
    hfid=fopen(histfile, 'a+');
end
%%
%define measurements to analyse based on if previous data has been found
if loadflag ==0
    looprange = 1:length(spec_name);
else
    if isequal(size(freq,1),size(ring_times,2),size(ring_amps,2),size(eignmode,2),size(Temp,1),size(phi,2),size(ia,2),ii)

        looprange = ii:length(spec_name);
    else %  program failed mid run - and some variables are not the right size.
        looprange = median([size(freq,1),size(ring_times,2),size(ring_amps,2),size(eignmode,2),size(Temp,1),size(phi,2),size(ia,2)]);
        % check for log file issues to avoid doubles 
        looprange =  looprange:length(spec_name);
        rawText = tex_import(histfile);
        if sum(~cellfun(@isempty,regexp(rawText,'scan_*','Match'))) == looprange(1)
            %delete last line of text file
            tmpfile = fullfile(FolderName,strcat('TMP_Log_',samplename,'_',datestr(now,'DD-mm-YYYY'),'.txt'));
            hfid2=fopen(tmpfile, 'a+');
            fprintf(hfid2,'Operator ID:\t %s' , ComputerID)
            fprintf(hfid2,'\nCreation Time:\t %s',datestr(now))
            fprintf(hfid2,'\nCurrent Directory:\t %s\n\n',pwd)
            fprintf(hfid2,'\nScanFileName\tRingfileName\tTemperature\tFrequency\tLabviewMode\tTau\tphi\tRejected\tRawRingdownLength\tOutputLength\txStart\txEnd\tGoodRingdown?\n')
            goodidx =find(~cellfun(@isempty,regexp(rawText,'scan_*','Match')));

            for k =1:numel(goodidx)-1
                fprintf(hfid2,rawText{goodidx(k)});
            end
            fclose(hfid2);
            [~,xfile] = fileparts(histfile);
            copyfile(histfile,fullfile(FolderName,strcat(xfile,'_old.txt')))
            copyfile(tmpfile,histfile)
            delete(tmpfile)
        end
        clear rawText hfid2 goodidx xfile k
        
        %check to see if variables are not the right size on reloading,
        %loop itteration will start again on median of variable sizes. 
        loopvars = {'freq','ring_times','ring_amps','eignmode','Temp','phi','ia'};
        varidx = find([size(freq,1),size(ring_times,2),size(ring_amps,2),size(eignmode,2),size(Temp,1),size(phi,2),size(ia,2)] ~=looprange(1));
        
        for kk = 1:numel(varidx)
            %check for right dimensions is Nx1 or 100000xN 
            if eval(sprintf('size(%s,1)>size(%s,2)',loopvars{varidx(kk)} ,loopvars{varidx(kk)} ))
                eval(sprintf('%s(end) = [];',loopvars{varidx(kk)}))
            else 
                 eval(sprintf('%s(:,end) = [];',loopvars{varidx(kk)}))
            end
        end 

        sanity = isequal(size(freq,1),size(ring_times,2),size(ring_amps,2),size(eignmode,2),size(Temp,1),size(phi,2),size(ia,2) );
        if ~sanity
            error('An error has been encountered trying to resume from current data')
            cprintf('err','\n Please seek assistance')
            return
        end 

    end
end
clear ii

for ii = looprange
    try
        %   check dimensions of tau
        if ii>1
            if size(tau,1)>1
                cprintf('err','\n Tau dimensions Changed:\tLoop Itt:%d\',ii)
                break
            end
        end

        try
            reload = ii;
            specdata = load(strcat(char(prefix),char(spec_name{ii})));
            ringdata = load(strcat('ring_',char(ring_name{ii})));
            ringdata2 = ringdata(:,2) ;
            ringdata1  =ringdata(:,1);
        catch
            error('Could not load in ringdown/scan_file data')
        end
        try
            freq(ii,:)  =  max(specdata(:,1));
        catch
            cprintf('err','\n Error: Scan file for %s is empty\n')
            freq(ii,:)  =  NaN;

        end
        SFilename = strcat(char(prefix),char(spec_name{ii}));
        RFilename = strcat('ring_',char(ring_name{ii}));

        if talk == 1
            fprintf('scan file being loaded: %s \n ',char(spec_name(ii,:)));
            fprintf('ring file being loaded: %s \n ',char(ring_name(ii,:)));
            fprintf('Mode frequency        : %.4f Hz \n',freq(ii,:));
        end
        %% extract mode numbers fron filenames
        test=0;
        while test==0
            if ii>1
                reload = ii;
            end

            ModeNumber = regexp(ring_name{ii},'\d{1,2}','Match');

            outputfile=strcat(FolderName,filesep,sprintf('Processed_Mode%d_',str2num(ModeNumber{oldnew})),samplename,'_',datestr(now,'DD-mm-YYYY'),'.txt');
            modes(ii) =str2num(ModeNumber{oldnew});


            %% extract temperatures from log files
            if ~exist('overtemp','var')

                try
                    % find log file for ringdown to extract temperature values
                    A = readtable(sprintf('log_m%d.txt',str2num(ModeNumber{oldnew})));
                    %find specific line in log file from ringdown file name
                    temp_idx = first(find(A.Var1 == str2num(ModeNumber{1})  & A.Var2 == str2num(ModeNumber{3})));
                    %extract temperature value
                    Temp(ii,:) = A.Var3(temp_idx);
                catch
                    Temp(ii,:)  = [];
                end
                if ~isempty(temp_idx)
                    Temp(ii,:) = A.Var3(first(find(A.Var1 == str2num(ModeNumber{1})  & A.Var2 == str2num(ModeNumber{3}))));
                else
                    Temp(ii,:) = NaN;
                    if talk
                        cprintf('err','No Matching Data for %s could be found in %s',ring_files1(ii).name,sprintf('log_m%d.txt',str2num(ModeNumber{oldnew})));
                        cprintf('err','\n Temperature : NaN - No Log File Entry Could be Found\n');
                    end
                    clear temp_idx
                end
            else
                Temp(ii,:) = str2num(cell2mat(overtemp));
            end

            %% perform intial fitting - prompt user
            if exist('auto','var') && ~isempty(auto)
                suppress = 'true';
                [tau(ii),phi(ii),gof(ii),reject(ii)]=curveAnalysis9(RFilename,SFilename,freq(ii,:),reload,suppress);

            else
                [tau(ii),phi(ii),gof(ii),reject(ii)]=curveAnalysis9(RFilename,SFilename,freq(ii,:),reload);
            end

            if exist('auto','var') && ~isempty(auto) && reject(ii ) == 0
                ShouldKeep= '' ;
                if  gof(ii)==0
                    break
                end
            elseif exist('auto','var') && ~isempty(auto) && reject(ii ) == 1
                ShouldKeep ='n';
                [tau(ii),phi(ii),gof(ii),reject(ii ),ringdata1,ringdata2] = curveAnalysis9(RFilename,SFilename,freq(ii,:),reload,'');
                if gof(ii) ==0
                    close
                end

            else
                suppress ='' ;
                if reject(ii) == 1
                    ShouldKeep = 'n';
                else
                    ShouldKeep = input('Leave empty to keep, 1 to re-fit, n to reject or r to re-do. \n','s');
                end
            end




            %% perform windowing if propmted
            if isempty(ShouldKeep)
                Keep(ii) = 1;
            elseif ShouldKeep=='1'
                [tau(ii),phi(ii),gof(ii),times,Amplitudes,goodpoints,cutpoints]=curveAnalysis3(RFilename,SFilename,freq(ii,:));
                newtimes = times(goodpoints>0);
                newAmplitudes = Amplitudes(goodpoints>0);

                Keep(ii) = 1;
                counter = counter +1 ;
            else
                Keep(ii) = 0;
            end


            fitagain='n';
            if isempty(fitagain)
                fitagain='n';
            end

            if ~exist('auto','var')
                close all
            end

            Q =exist('newtimes');
            if Q
                if exist('auto','var') && ~isempty(auto)
                    eignmode(:,ii)   = str2num(ModeNumber{oldnew});
                else
                    % make sure that each array starts at index zero
                    if newtimes(1)>0
                        newtimes = newtimes-newtimes(1);
                    end
                    ring_times(:,ii) = [newtimes; NaN(matchk-length(newtimes),1)];
                    ring_amps(:,ii)  = [newAmplitudes; NaN(matchk-length(newAmplitudes),1)];
                    ia(:,ii)         = [goodpoints(goodpoints>0) ;NaN(matchk-length(goodpoints(goodpoints>0)),1)];
                    eignmode(:,ii)   = str2num(ModeNumber{oldnew});
                end

                fid=fopen(outputfile, 'a+');
                fprintf(fid,'%2.4f\t%5.4f\t%2.3f\t%g\t%d\n',...
                    Temp(ii),freq(ii),tau(ii),gof(ii),Keep(ii));
                fclose(fid);

                if floor(cutpoints(1))<0
                    cutpoints(1) = 0 ;
                end


                if talk == 1
                    fprintf('*Manually -Refit*\n');
                    fprintf('length of original file  Ring: %d ,  spec: %d\n',length(ringdata2),length(ringdata1));
                    fprintf('length of new file       Ring: %d ,  spec: %d\n',length(newAmplitudes),length(newtimes));
                    fprintf('ones: %d, zeros, %d, NaNs %d \n',numel(find(ia(:,ii))),numel(find(ia(:,ii)==0)),numel(find(ia(:,ii)==NaN)));
                    fprintf('\n')
                end
                test=1 ;
            end

            if Keep(ii) == 0
                if ShouldKeep=='n'
                    fprintf('*Ringdown rejected*\n')

                    if exist('auto','var') && ~isempty(auto)
                        eignmode(:,ii)   = str2num(ModeNumber{oldnew});
                        freq(ii,:)       = freq(ii,:);

                    else
                        ring_times(:,ii) = [NaN(length(ringdata1),1); NaN(matchk-length(ringdata1),1)];
                        ring_amps(:,ii)  = [NaN(length(ringdata2),1); NaN(matchk-length(ringdata2),1)];
                        tmp =[NaN(length(ringdata2),1); NaN(matchk-length(ringdata2),1)];
                        eignmode(:,ii)   = NaN;
                        freq(ii)       = NaN;
                        ia(:,ii)        = tmp ;
                        tau(ii)        = NaN;
                        phi(ii)        = NaN;
                        Temp(ii)       = NaN;

                    end

                    fid=fopen(outputfile, 'a+');
                    fprintf(fid,'%2.4f\t%5.4f\t%2.3f\t%g\t%d\n',...
                        Temp(ii),freq(ii),tau(ii),gof(ii),Keep(ii));
                    fclose(fid);

                    if talk == 1
                        fprintf('length of original file  Ring: %d ,  spec: %d\n',length(ringdata2),length(ringdata1));
                        fprintf('length of new file       Ring: %d ,  spec: %d\n',length(ringdata2 ),length(ringdata1));
                        fprintf('Ringdown rejected - not written to Python output' );
                        fprintf('\n');

                    end
                    test=1;
                end
            elseif Keep(ii) == 1 && test ==0
                if exist('auto','var') && ~isempty(auto)
                    eignmode(:,ii)   = str2num(ModeNumber{oldnew}) ;

                else
                    ring_times(:,ii) = [ringdata1; NaN(matchk-length(ringdata1),1)];
                    ring_amps(:,ii)  = [ringdata2; NaN(matchk-length(ringdata2),1)];
                    eignmode(:,ii)   = str2num(ModeNumber{oldnew}) ;

                    if exist('goodpoints','var')
                        ia(:,ii)         = [goodpoints(goodpoints>0); NaN(matchk-length(goodpoints(goodpoints>0),1))];
                    else
                        ia(:,ii)         = [ones(size(ringdata2)); NaN(matchk-length(ones(size(ringdata2))),1)];
                    end

                end

                fid=fopen(outputfile, 'a+');
                fprintf(fid,'%2.4f\t%5.4f\t%2.3f\t%g\t%d\n',...
                    Temp(ii),freq(ii),tau(ii),gof(ii),Keep(ii));
                fclose(fid);

                if talk == 1

                    if exist('auto','var') && ~isempty(auto)
                        fprintf('length of original file  Ring: %d ,  spec: %d\n',length(ringdata2),length(ringdata1));
                        fprintf('length of new file       Ring: %d ,  spec: %d\n',length(ringdata2),length(ringdata1));
                        fprintf('\n');
                    else
                        fprintf('length of original file  Ring: %d ,  spec: %d\n',length(ringdata2),length(ringdata1));
                        fprintf('length of new file       Ring: %d ,  spec: %d\n',length(ringdata2),length(ringdata1));
                        fprintf('ones: %d, zeros, %d, NaNs %d \n',numel(find(ia(:,ii))),numel(find(ia(:,ii)==0)),numel(find(ia(:,ii)==NaN)));
                        fprintf('\n');
                    end
                end
                test=1;
            elseif strcmp(fitagain,'y')
                C=C+1;
            end

            if ShouldKeep == 'r'
                if ii >1
                    ii = reload-1;
                    clear reload
                end
            end
            if ~isequal(size(freq,1),size(ring_times,2),size(ring_amps,2),size(eignmode,2),size(Temp,1),size(phi,2),size(ia,2))
                cprintf('err','\n Internal Error - Size Mismatch on outputs detected\n')
                cprintf('err','\n Loop Itteration :%d',ii)
                cprintf('err','\n Ringdown File :%s',strcat(char(prefix),char(spec_name{ii})))
                cprintf('err','\n Scan File: %s',strcat('ring_',char(ring_name{ii})))
                cprintf('err','\n Saving Current state...\n')
                save(fullfile(FolderName,'Recovered_outputs.mat'))
                cprintf('err','\n Exciting...\n')
                return
            end

            %write outputs to logfile
            if ii ==1
                fprintf(hfid,'Operator ID:\t %s' , ComputerID)
                fprintf(hfid,'\nCreation Time:\t %s',datestr(now))
                fprintf(hfid,'\nCurrent Directory:\t %s\n\n',pwd)
                fprintf(hfid,'\nScanFileName\tRingfileName\tTemperature\tFrequency\tLabviewMode\tTau\tphi\tRejected\tRawRingdownLength\tOutputLength\txStart\txEnd\tGoodRingdown?')
            end

            if exist('cutpoints','var')
                a = table(strcat(prefix,char(ring_name{ii})),{strcat('ring_',char(ring_name{ii}))},Temp(ii),freq(ii),modes(ii),tau(ii),phi(ii),reject(ii),length(newAmplitudes),length(ringdata1),floor(cutpoints(1)),floor(cutpoints(2)),Keep(ii));
                writetable(a,histfile,'Delimiter','\t','WriteMode','append')
            else
                a = table(strcat(prefix,char(ring_name{ii})),{strcat('ring_',char(ring_name{ii}))},Temp(ii),freq(ii),modes(ii),tau(ii),phi(ii),reject(ii),length(ringdata1),length(ringdata1),NaN,NaN,Keep(ii));
                writetable(a,histfile,'Delimiter','\t','WriteMode','append')
            end

            clear goodpoints ringdata ringdata1 ringdata2 specdata Q tmp newtimes newAmplitudes ShouldKeep cutpoints a


        end

        if isnan(eignmode(ii)) && strcmpi(Keep(ii),'n')
            cprintf('err','\nItteration %d Gives NAN eignmode\n',ii)
            dbstop
        end



        clear specdata ringdata ringdata2 fitresult b a cutpoints
    catch
        save(fullfile(FolderName,'Recovered_outputs.mat'))
        ii = max(looprange);

    end
end


%pythondata
pyname = strcat(FolderName,'/PyTotalAnalysis.mat');

fprintf('\n Analysis Summary:\n**************************\n')

fprintf('Valid Ringdowns         :\t%d\n',numel(find(Keep )))
fprintf('Ringdowns Manually Refit:\t%d\n', counter)
fprintf('Rejected Ringdowns      :\t%d\n\n', numel(find(Keep ==0)))

rough_f =unique(sort(ceil(freq(Keep'>0 & ~isnan(freq))./100)*100));
counts = groupcounts(sort(ceil(freq(Keep'>0 & ~isnan(freq))./100)*100));

fprintf('Rough Freq\t\tNo.Measurements\n')

for k = 1:numel(rough_f)
    if rough_f(k)<10000
        fprintf('  %2.2f%+22s\n',rough_f(k),string(counts(k)))
    else
        fprintf('  %2.2f%+21s\n',rough_f(k),string(counts(k)))
    end
end


if ~exist('auto','var')
    freq =freq(find(~isnan(freq)&freq>0));
    modes= eignmode(find(~isnan(freq)&freq>0))';
    tau = [tau(find(~isnan(freq(:,1))&freq(:,1)>0))];
    phi = [phi(find(~isnan(freq(:,1))&freq(:,1)>0))];
    Temp= [Temp(find(~isnan(freq)&freq>0))];
    ring_amps =ring_amps(:,any(ring_amps));
    ring_times =ring_times(:,any(ring_times));
    ia = ia(:,find(sum(~isnan(ia))>0));

    save(pyname,'freq','ring_times','ring_amps','modes','ia','tau','phi','Temp','-v7.3');

else
    freq =freq(find(~isnan(freq)&freq>0));
    modes= eignmode(find(~isnan(freq)&freq>0))';
    tau = [tau(find(~isnan(freq(:,1))&freq(:,1)>0))];
    phi = [phi(find(~isnan(freq(:,1))&freq(:,1)>0))];
    Temp= [Temp(find(~isnan(freq(:,1))&freq(:,1)>0))];
    freq = freq(~isnan(freq))
    reject =~reject;
    save(pyname,'freq','modes','tau','phi','gof','reject','Temp','-v7.3');
    data = table(Temp,freq,tau', phi' , gof', reject', eignmode','VariableNames',{'Temp','Freq','tau','Loss','GOF','Keep','Mode',});

    colors = distinguishable_colors(max(eignmode));
    figure
    hold on
    grid on
    for i =1: max(eignmode)
        eval(sprintf('mode%d = data(find(eignmode ==i),:);',i))
        if  ~isempty(eval(sprintf('mode%d',i)))
            writetable(eval(sprintf('mode%d(:,1:end-1)',i)), sprintf('Processed_Mode_%d.txt',i),'Delimiter','\t','WriteVariableNames',0);

            eval(sprintf('mode%d = sortrows(mode%d,1);',i,i))
            eval(sprintf('plot(mode%d.Temp,mode%d.Loss,''o'',''DisplayName'',''Mode %d'')',i,i,i))

        end
    end
    title(regexprep(samplename,'_',' '))
    set(gca,'YScale','log')
    legend('Location','BestOutside')
    xlabel('Temperature [K]')
    ylabel('Coated Loss [ \phi_{coating} ]')
    savefig(gcf,strcat(samplename,'.fig'))

end

if isfile(fullfile(FolderName,'Recovered_outputs.mat'))
    delete(fullfile(FolderName,'Recovered_outputs.mat'))
end

%run python analysis
cd(FolderName);

%functionallity currently uses python and Linux base commands - this will
%not work for other systems and needs to be fixed before rolling out to
%everyone
if contains(ComputerID,'ComputerID')
    reply = input('\nPerform Python Analysis? [Y/N] \n','s');
    if contains(lower(reply),'y')

        pre_fit_py

        pause(2)
        system(sprintf('! open %s ',char(fnames('results*.png'))))

        reply2 = input('\n Do you want to go through this data now?  [Y/N]','s');

        if contains(lower(reply2),'y')|| isempty(reply2)
            beat_excel
        end
    end
end
end

%% Checking for new versions and self-updating
% Function used to fetch latest master-branch from github and install in
% MATLAB path. Code taken adapted from export_fig.m https://github.com/altmany/export_fig/

% Update the installed version of Prebeat from the latest version online
function updateInstalledVersion()
% Download the latest version of Prebeat into the Prebeat folder
zipFileName = 'https://github.com/Titian2/Prebeat/archive/master.zip';
fprintf('Downloading latest version of %s from %s...\n', 'Prebeat', zipFileName);
folderName = fileparts(which(mfilename('fullpath')));
targetFileName = fullfile(folderName, datestr(now,'yyyy-mm-dd.zip'));
try
    folder = hyperlink(['matlab:winopen(''' folderName ''')'], folderName);
catch  % hyperlink.m is not properly installed
    folder = folderName;
end
try
    urlwrite(zipFileName,targetFileName); %#ok<URLWR>
catch err
    error('Prebeat_update:download','Error downloading %s into %s: %s\n',zipFileName,targetFileName,err.message);
end

% Unzip the downloaded zip file in the Prebeat folder
fprintf('Extracting %s...\n', targetFileName);
try
    unzip(targetFileName,folderName);
    % Fix issue #302 - zip file uses an internal folder Prebeat-master
    subFolder = fullfile(folderName,'Prebeat-master');
    try movefile(fullfile(subFolder,'*.*'),folderName, 'f'); catch, end %All OSes
    try movefile(fullfile(subFolder,'*'),  folderName, 'f'); catch, end %MacOS/Unix
    try movefile(fullfile(subFolder,'.*'), folderName, 'f'); catch, end %MacOS/Unix
    try rmdir(subFolder); catch, end
catch err
    error('Prebeat:update:unzip','Error unzipping %s: %s\n',targetFileName,err.message);
end

% Notify the user and rehash
fprintf('Successfully installed the latest %s version in %s\n', mfilename, folder);
clear functions %#ok<CLFUNC>
rehash
end


% Check for newer version (only once a day)
function isNewerVersionAvailable = checkForNewerVersion_MATLAB(currentVersion)
persistent lastCheckTime lastVersion
isNewerVersionAvailable = false;
%     if nargin < 1 || isempty(lastCheckTime) || now - lastCheckTime > 1
url ='https://raw.githubusercontent.com/Titian2/Prebeat/main/prebeat.m';
try
    str = webread(url);
    cll = strsplit(str, '\n')';
    versionHistory = regexp(str, '\%\s+v\d+[.](.*?)\n','Match')';
    versionHistory = regexp(versionHistory,'(v\d\.\d) | (v\d\.\d\w) | (v\d\.\d\.\d)','Match');
    versionHistory = [versionHistory{:}]';
    versionHistory = cellfun(@(str)regexprep(str,'\s+', ''),versionHistory ,'UniformOutput',false);

    a = cellfun(@(str)regexp(str,'\d+','Match'),versionHistory,'UniformOutput',false);

    for i =1:numel(a)
        if size(a{i},2) < 3
            a{i} = [a{i},'0'];
        end
        versionNum(i) = str2num(cell2mat(a{i}(1:3)));
    end

    latestVersion = max(versionNum)/100;

    if nargin < 1
        currentVersion = lastVersion;
    else
        currentVersion = currentVersion + 1e3*eps;
    end
    isNewerVersionAvailable = latestVersion > currentVersion;
    if isNewerVersionAvailable
        try
            versionDesc = cll(find(contains(cll,versionHistory)));
            versionDesc = regexprep(extractAfter(versionDesc,versionHistory),'\s{2,}','');
        catch
            % Something bad happened - only display the latest version description
            versionDesc = {'version description not avalible'};
        end
        try versionDesc = strjoin(strrep(strcat(' ***', strtrim(strsplit(versionDesc,';'))),'***','* '), char(10)); catch, end %#ok<CHARTEN>
        msg = sprintf(['You are using version %g of Prebeat. A newer version (%g) is available, with the following improvements/fixes:\n\n %s\n\n',... ,
            'A change-log of recent releases is available here; the complete change-log is included at the top of the Prebeat.m file.\n',...
            'You can download the new version from GitHub'],currentVersion, latestVersion,string(versionDesc{1}));
        msg = hyperlink('https://github.com/Titian2/Prebeat/', 'GitHub', msg);
        msg = hyperlink('matlab:Prebeat(''-update'')', msg);
        warning('Prebeat:version',msg);
    end
catch
    % ignore
end
lastCheckTime = now;
lastVersion = currentVersion;
end




