function [tau,phi,gof,reject,time,Amplitude]=curveAnalysis9(Rfilename,Sfilename,freq,reload,suppress)
%
%Function takes in ringdown,scan and an equivilent frequency value and
%produces an exponential fit to the input data. Calculating the mechanical
%loss from the exponential fitting
%
% [tau,phi,gof,reject]=curveAnalysis9(Rfilename,Sfilename,freq,reload)
%
% INPUTS
% Rfilename : filename of ringdown file as string
% Sfilename : filename of matching scan file as string
% freq      : frequency of the ringdown - used for Mechanical loss calc.
%
% OUTPUTS
% tau       : decay constant of the fitted ringdown
% phi       : mechanical loss calculated from tau
% gof       : details on the 'goodness of fit'
% reject    : binary value which dictates if measurement is to be rejected
%             from further analysis
%
%
% requires cbblindplot for plotting - plots in colorblind friendly colors
% v1.2.1 - Small changes to fix figures not being generated 05/12/2022
% v1.2   - Changing handling of 'Not Enough Data' error
%          :S.Tait 03/2021
%
%
% - Edited by S.Tait 2022
% - s.tait.1@research.gla.ac.uk
%


currentVersion = 1.2;
noscan = 0;

cb=cblindplot;
scrsz = get(0,'ScreenSize');
XXX = [1 1*scrsz(4)/6 scrsz(3) 2*scrsz(4)/3];
%%%%load data
A=load(Rfilename);
B=load(Sfilename);

% RLength=size(A,1);
% FivePercent=ceil(0.05*RLength)+1;
%
time=A(:,1);
Amplitude=A(:,2);


if ~isempty(B)
    time2=B(:,1);
    Amplitude2=B(:,2);
else
    noscan = 1;
end



% %clean spikes from ringdown data
time(find((Amplitude>2000))) = NaN;
Amplitude(find((Amplitude>2000))) = NaN;




%%%fit data
if length(time)<=3
    tau=0;
    phi=0;
    gof=0;
    
    figure;
    plot(1:5,1:5,1:5,5:-1:1);
    text(3,3,'Not Enough Data','HorizontalAlignment','Center');
    reject = 1;
else
    
    [D,E]=fit(time(~isnan(time)),Amplitude(~isnan(Amplitude)),'exp1');
    
    if D.b>0
        reject = 1 ;
    else 
        reject = 0; 
    end

    if ~exist('suppress','var')
        if  reject ==0
            figure('Position',XXX);
            subplot(1,2,1);
            hold on
            grid on
            plot(time,Amplitude,'bo')
            plot(time,D(time),'color',cb(2,:),'LineWidth',4);
            
            xlabel('Time ')
            ylabel('Amplitude')
            if numel(isnan(time))>1
                legend('Cleaned : Spikes Removed ')
            end
            
            set(subplot(1,2,1),'FontSize',14)
            
            subplot(1,2,2);
            
            hold on
            grid on
            if noscan == 0
                plot(time2,Amplitude2,'color',cb(4,:),'LineWidth',3);
                
                if numel(isnan(time2))>1
                    legend('Cleaned : Spikes Removed ')
                end
                
            else
                plot(1:5,1:5,1:5,5:-1:1);
                text(3,3,'Not Enough Data','HorizontalAlignment','Center');
            end
            
            
            xlabel('Frequency ')
            ylabel('Amplitude')
            set(subplot(1,2,2),'FontSize',14)           
            title(sprintf('%0.3f Hz',freq),'FontSize',20)
            pause(0.5) % for some reason systems seem to  slow down
            % or not plot the figures at all without this pause in place
        end
    end
    
    %%%report data
    tau=-1/D.b;
    
    if isnan(tau) 
        cprintf('err','\nERR: tau is NAN\n') 
        pause 
    end 
        
    phi=1/(pi*freq*tau);
    gof=E.adjrsquare;
    
end

if isempty(reject)
    disp('stop')
end

if sum(isnan(time))>0
    time =time - first(time(~isnan(time)));
end



end