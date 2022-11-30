function [tau,phi,gof,time,Amplitude,ia,x]=curveAnalysis3(Rfilename,Sfilename,freq)
%
%Function takes in ringdown,scan and an equivilent frequency value and
%produces an exponential fit to the input data. Calculating the mechanical
%loss from the exponential fitting. Allows for windowing and refitting of
%data 
%
%  [tau,phi,gof,newtime,newAmplitudes,ia]=curveAnalysis3(Rfilename,Sfilename,freq)
%
% INPUTS 
% Rfilename     : filename of ringdown file as string 
% Sfilename     : filename of matching scan file as string 
% freq          : frequency of the ringdown - used for Mechanical loss calc. 
%
% OUTPUTS 
% tau           : decay constant of the fitted ringdown 
% phi           : mechanical loss calculated from tau 
% gof           : details on the 'goodness of fit' 
% newtime       : vector of new time values ( if data has been re-fitted ) 
% newAmplitudes : vector of new Amplitude values ( if data has been re-fitted ) 
% ia            : boolean vector used for further python analysis 
%
% requires cbblindplot for plotting - plots in colorblind friendly colors 
% 
% - Edited by S.Tait 
% - s.tait.1@research.gla.ac.uk 
% 

cb=cblindplot;
scrsz = get(0,'ScreenSize');
XXX = [1 1*scrsz(4)/6 scrsz(3) 2*scrsz(4)/3];
%%%%load data
A=load(Rfilename);
B=load(Sfilename);

time=A(:,1);
Amplitude=A(:,2);

time2=B(:,1);
Amplitude2=B(:,2);

f1=figure('Position',XXX);
hold on;
plot(time,Amplitude,'.');

[D,E]=fit(time,Amplitude,'exp1');
fity=D.a*exp(D.b*time);

plot(time,Amplitude,'b.',time,fity,'r-');  
E.adjrsquare
%%%check data with user
test=0;

while(test==0)
      
    
    %if strcmp(B,'n')
    gca ; 
    title('Select region to use')
    [x,y]=ginput(2);
     
    xstart=x(1);
    xstop=x(2);
        
    C=1;
    for n=1:length(time)
        if (time(n)>=xstart) && (time(n)<=xstop)
            newtime(C)=time(n);
            newAmplitudes(C)=Amplitude(n);
            C=C+1;
        end
    end
    newAmplitudes = newAmplitudes';
    newtime = newtime';
    
    f1=figure('Position',XXX);
    hold on;
    plot(newtime,newAmplitudes,'.');
    [D,E]=fit(newtime,newAmplitudes,'exp1');
    fity=D.a*exp(D.b*newtime);
    
    plot(newtime,newAmplitudes,'b.',newtime,fity,'r-');  
    E.adjrsquare
    B=input('Keep this fit?  y/n/q [y] \n','s');

    if isempty(B)
        B='y';
    end
    if strcmp(B,'y')
        test=1;
    elseif strcmp(B,'q')
        time = 0;
        test = 1;
    elseif strcmp(B,'n')
        test = 0;
    else
        fprintf('please, only y, n or q (lower case).  \n');
    end
    close(f1);
end

%%%fit data
if length(time)<=3
    tau=0;
    phi=0;
    gof=0;
    
    figure;
    plot(1:5,1:5,1:5,5:-1:1);
    text(3,3,'Not Enough Data','HorizontalAlignment','Center');
else

    [D,E]=fit(newtime,newAmplitudes,'exp1');

    fity=D.a*exp(D.b*newtime);

    figure('Position',XXX);
    subplot(1,2,1);plot(newtime,newAmplitudes,'b.',newtime,fity,'color',cb(1,:));
    subplot(1,2,2);plot(time2,Amplitude2,'color',cb(4,:));
    xlabel('Time ') 
    ylabel('Amplitude')
    %return qualifer for python script
    ia = abs(double(time<x(1)|time>x(2))-1);

%%%report data
    tau=-1/D.b;
    phi=1/(pi*freq*tau);
    gof=E.adjrsquare;
    
  
    
    
end