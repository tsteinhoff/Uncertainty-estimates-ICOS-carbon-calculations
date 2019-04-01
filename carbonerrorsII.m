%%
tic
close all

dataout=[];

FCO2=[250, 300, 350, 400, 450, 500, 550, 600, 700, 800];   % i parameter
T=[5, 10, 15, 20, 25];                           % j parameter

SAL=35;
PRESIN=5; PRESOUT=5; SI=0; PO4=0; 
pHSCALEIN=1;
K1K2CONSTANTS=10;
KSO4CONSTANTS=1;
eSAL= 0.1;
eTEMP= 0.005;
eSI= 0;
ePO4= 0;
epK=[0.002, 0.0055, 0.01, 0.01, 0.01, 0.02, 0.02];
eBt= 0.02;
r= 0;

efco2=2; %0.5;
ddic=2; %[1.5 3 4.5];
dta=2; %[2, 4, 6];
dph=0.001;

eDIC=ddic;
eTA=dta;
epH=dph;

%

for i=1:length(FCO2)
    % establishing consitent carbonsystem @ 20°C with DIC=1950 µmol/kg
    fco2i=FCO2(i);
    DIC=1950;
           
    for j=1:length(T)
        temp=T(j);
        % fCO2 at temperature of interest
        %[RESULT,HEADERS,NICEHEADERS]=CO2SYS(fco2_20,DIC,5,2,SAL,20,temp,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        %fco2_temp=RESULT(20);
        
        % carbonate system at temperature of interest
        [RESULT,HEADERS,NICEHEADERS]=CO2SYS(fco2i,DIC,5,2,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        fco2=fco2i;
        dic=DIC;
        ta=RESULT(1);
        ph=RESULT(3); H=10^(-ph);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calc TA and pH as f(fCO2, DIC)
        PAR1=fco2; PAR2=DIC; PAR1TYPE=5; PAR2TYPE=2; ePAR1=efco2; ePAR2=eDIC;
        [err, ehead, eunits] =      errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        dta_fco2_dic=err(1);
        eh=err(2)*1e-9;  
        dph_fco2_dic=-log10(1+eh/H);
        
        % calc pH and DIC as f(fCO2, TA)
        PAR1=fco2; PAR2=ta; PAR1TYPE=5; PAR2TYPE=1; ePAR1=efco2; ePAR2=eTA;
        [err, ehead, eunits] =      errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        ddic_fco2_ta=err(1);
        eh=err(2)*1e-9; 
        dph_fco2_ta=-log10(1+eh/H);
        
        % calc TA and DIC as f(fCO2, pH)
        PAR1=fco2; PAR2=ph; PAR1TYPE=5; PAR2TYPE=3; ePAR1=efco2; ePAR2=epH;
        [err, ehead, eunits] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        ddic_fco2_ph=err(2);
        dta_fco2_ph=err(1);
        
        % calc fCO2 and pH as f(TA, DIC)
        PAR1=ta; PAR2=dic; PAR1TYPE=1; PAR2TYPE=2; ePAR1=eTA; ePAR2=eDIC;
        [err, ehead, eunits] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        %fprintf('%s   %s %s %s %s  %s %s %s %s \n', ehead{1:9});
        %fprintf('%s   %s  %s    %s %s %s    %s      %s          %s \n', eunits{1:9});
        %fprintf('%f  %f  %f  %f  %f %f  %f     %f     %f \n', err(1:9));
        dfco2_ta_dic=err(3);
        eh=err(1)*1e-9; 
        dph_ta_dic=-log10(1+eh/H);
        
        % calc fCO2 and TA as f(pH, DIC)
        PAR1=ph; PAR2=dic; PAR1TYPE=3; PAR2TYPE=2; ePAR1=epH; ePAR2=eDIC;
        [err, ehead, eunits] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
%         fprintf('%s   %s %s %s %s  %s %s %s %s \n', ehead{1:9});
%         fprintf('%s   %s  %s    %s %s %s    %s      %s          %s \n', eunits{1:9});
%         fprintf('%f  %f  %f  %f  %f %f  %f     %f     %f \n', err(1:9));
        dfco2_pH_dic=err(3);
        dta_pH_dic=err(1);
        
        % calc fCO2 and DIC as f(pH, TA)
        PAR1=ph; PAR2=ta; PAR1TYPE=3; PAR2TYPE=1; ePAR1=epH; ePAR2=eTA;
        [err, ehead, eunits] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,temp,temp,PRESIN,PRESOUT,SI,PO4,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
%         fprintf('%s   %s %s %s %s  %s %s %s %s \n', ehead{1:9});
%         fprintf('%s   %s  %s    %s %s %s    %s      %s          %s \n', eunits{1:9});
%         fprintf('%f  %f  %f  %f  %f %f  %f     %f     %f \n', err(1:9));
        dfco2_pH_ta=err(3);
        ddic_pH_ta=err(1);
        
        
        
        
        
        
        
        
        
        
        
        dataout=[dataout; temp, fco2, dta_fco2_dic, dph_fco2_dic, ddic_fco2_ta, dph_fco2_ta, ddic_fco2_ph, dta_fco2_ph];
   
    end
end
%

temp=dataout(:,1);
fco2=dataout(:,2);
dta_fd=dataout(:,3);
dph_fd=dataout(:,4);
ddic_ft=dataout(:,5);
dph_ft=dataout(:,6);
ddic_fp=dataout(:,7);
dta_fp=dataout(:,8);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
minfco2=200;
maxfco2=800;
toc
figure
close all
%%

close all
figure
x0=10;
y0=10;
width=1300;
height=650;
set(gcf,'position',[x0,y0,width,height])
% plots f(fCO2, pH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TA
subplot(3,2,1)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),dta_fp(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
%xlabel('fCO_2 [µatm]')
ylabel('u(TA) [µmol/kg]')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['TA calculated from fCO_2 and pH'];
title(ti)

y1 = 2;
line([minfco2,maxfco2],[y1,y1],'Color','r')
y2=10;
line([minfco2,maxfco2],[y2,y2])
text(minfco2+5,y1+2,'Climate goal','Color','r')
text(minfco2+5,y2+2,'Weather goal','Color','b')
%text(minfco2+5,y1+5,'Climate goal','Color','r')
%text(minfco2+455,y2+5,'Weather goal','Color','b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot DIC
subplot(3,2,3)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),ddic_fp(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
%xlabel('fCO_2 [µatm]')
ylabel('u(DIC) [µmol/kg]')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['DIC calculated from fCO_2 and pH'];
title(ti)

y = 2;
line([minfco2,maxfco2],[y,y],'Color','r')
y=10;
line([minfco2,maxfco2],[y,y])

% plots f(fco2,DIC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TA
%figure(2)
subplot(3,2,2)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),dta_fd(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
%xlabel('fCO_2 [µatm]')
ylabel('u(TA) [µmol/kg]')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['TA calculated from fCO_2 and DIC'];
title(ti)

y = 2;
line([minfco2,maxfco2],[y,y],'Color','r')
y=10;
line([minfco2,maxfco2],[y,y])
ylim([0, 12])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot pH
subplot(3,2,5)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),-dph_fd(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
  
xlabel('fCO_2 [µatm]')
ylabel('u(pH)')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['pH calculated from fCO_2 and DIC'];
title(ti)


y = 0.003;
line([minfco2,maxfco2],[y,y],'Color','r')
y=0.02;
line([minfco2,maxfco2],[y,y])
ylim([0, 0.025])
%   yt = get(gca,'YTick');
%     set(gca,'YTickLabel', sprintf('%.4f\n',yt));

    
    
    
% plot f(fco2, TA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot DIC
subplot(3,2,4)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),ddic_ft(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
%xlabel('fCO_2 [µatm]')
ylabel('u(DIC) [µmol/kg]')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['DIC calculated from fCO_2 and TA'];
title(ti)

y = 2;
line([minfco2,maxfco2],[y,y],'Color','r')
y=10;
line([minfco2,maxfco2],[y,y])
ylim([0, 12])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot pH
subplot(3,2,6)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),-dph_ft(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');

xlabel('fCO_2 [µatm]')
ylabel('u(pH)')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['pH calculated from fCO_2 and TA'];
title(ti)
% % 
y = 0.003;
line([minfco2,maxfco2],[y,y],'Color','r')
y=0.02;
line([minfco2,maxfco2],[y,y])
ylim([0, 0.025])
%   yt = get(gca,'YTick');
%   set(gca,'YTickLabel', sprintf('%.3f\n',yt));







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% plot only f(DIC)
% plot TA from DIC
subplot(1,2,1)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),dta_fd(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
xlabel('fCO_2 [µatm]')
ylabel('u(TA) [µmol/kg]')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['TA calculated from fCO_2 and DIC'];
title(ti)
% plot pH from DIC
subplot(1,2,2)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),-dph_fd(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
  yt = get(gca,'YTick');
  set(gca,'YTickLabel', sprintf('%.4f\n',yt));
xlabel('fCO_2 [µatm]')
ylabel('u(pH)')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['pH calculated from fCO_2 and DIC'];
title(ti)









%% plot only f(TA)
% plot DIC
figure
subplot(1,2,1)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),ddic_ft(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
xlabel('fCO_2 [µatm]')
ylabel('u(DIC) [µmol/kg]')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['DIC calculated from fCO_2 and TA'];
title(ti)

% plot pH
subplot(1,2,2)
scatter(fco2(fco2>minfco2 & fco2<maxfco2),-dph_ft(fco2>minfco2 & fco2<maxfco2),20,temp(fco2>minfco2 & fco2<maxfco2),'filled');
 yt = get(gca,'YTick');
 set(gca,'YTickLabel', sprintf('%.4f\n',yt));
xlabel('fCO_2 [µatm]')
ylabel('u(pH)')
grid on
h = colorbar;
set(get(h,'title'),'string','Temp [°C]');
ti=['pH calculated from fCO_2 and TA'];
title(ti)


