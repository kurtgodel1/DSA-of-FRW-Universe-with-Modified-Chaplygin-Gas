close



Fcolor=[0.5 0 1];
dScolor=.6*[0 0.7 0.5];
Mcolor=[0.7 0.2 0.2];
Ecolor=[1 0 0];


fontsize=20;
textboxfontsize=20; 
linewidth=2.5;

mx=0.5;
my=0.9;
fx=0.7;
fy=0.5;
dsx=1-fx;
dsy=fy;
ex=mx;
ey=0.1;

%lines
annotation('line',[mx fx],[my fy]);
annotation('line',[fx dsx],[fy dsy]);
annotation('line',[dsx mx],[dsy my]);
annotation('line',[fx ex],[fy ey]);
annotation('line',[ex dsx],[ey dsy]);

%texts


%M and around
txtM=annotation('textbox',[mx-0.05 my-0.03 0.1 0.1],'String','\textbf{M}','EdgeColor','none','HorizontalAlignment',"center",'interpreter','latex');
txtM.FontSize=fontsize;
txtM.Color=Mcolor;

textMf=annotation('textbox',[mx+0.07 my-0.11 0.1 0.1],'String',{'$H < 0$','$\Omega_\Lambda+\Omega_A = 0$'},'interpreter','latex');
textMf.EdgeColor=Mcolor;
textMf.FontSize=textboxfontsize;
textMf.LineWidth=linewidth;

textMCD=annotation('textbox',[mx-0.2 my-0.11 0.1 0.1],'String',{'$H > 0$','$\Omega_A = \sqrt[\alpha + 1]{\gamma}\Omega$'},'interpreter','latex');
textMCD.EdgeColor=Mcolor;
textMCD.FontSize=textboxfontsize;
textMCD.LineWidth=linewidth;


%F and around
txtF=annotation('textbox',[fx+0.01 fy-0.07 0.1 0.1],'String','\textbf{F}','EdgeColor','none','interpreter','latex');
txtF.FontSize=fontsize;
txtF.Color=Fcolor;

textFm=annotation('textbox',[fx-0.03 fy+0.08 0.1 0.1],'String',{'$H > 0$','$\Omega_\Lambda+\Omega_A = 0$'},'interpreter','latex');
textFm.EdgeColor=Fcolor;
textFm.FontSize=textboxfontsize;
textFm.LineWidth=linewidth;

textFCD=annotation('textbox',[fx-0.255 fy-0.1 0.1 0.1],'String',{'$H > 0$', 'k = 0','$\Omega_\Lambda+\Omega_A \neq 0$ if k = -1','$\Omega_\Lambda+\Omega_A > (\Omega_\Lambda+\Omega_A)_c$ if k = +1'},'interpreter','latex');
textFCD.EdgeColor=Fcolor;
textFCD.FontSize=textboxfontsize*0.6;
textFCD.LineWidth=linewidth;

textFf=annotation('textbox',[fx+0.13 fy-0.06 0.1 0.1],'String',{'$H > 0$', 'k = +1','$\Omega_\Lambda + \Omega_A < (\Omega_\Lambda + \Omega_A)_c$'},'interpreter','latex');
textFf.EdgeColor=Fcolor;
textFf.FontSize=textboxfontsize*.7;
textFf.LineWidth=linewidth;

textFe=annotation('textbox',[fx-0.04 fy-0.2 0.1 0.1],'String',{'$H > 0$','k=+1' '$\Omega_\Lambda + \Omega_A = (\Omega_\Lambda + \Omega_A)_c$'},'interpreter','latex');
textFe.EdgeColor=Fcolor;
textFe.FontSize=textboxfontsize*.7;
textFe.LineWidth=linewidth;


%dS and around
txtCD=annotation('textbox',[dsx-0.04 dsy-0.07 0.1 0.1],'String','\textbf{CD}','EdgeColor','none','interpreter','latex');
txtCD.FontSize=fontsize;
txtCD.Color=dScolor;

textCDf=annotation('textbox',[dsx+0.08 dsy+0.025 0.1 0.1],'String',{'$H < 0$', 'k = 0','$\Omega_A \neq \sqrt[\alpha + 1]{\gamma}\Omega$ if k = -1','$\Omega_\Lambda+\Omega_A < (\Omega_\Lambda+\Omega_A)_c$ if k = +1'},'interpreter','latex');
textCDf.EdgeColor=dScolor;
textCDf.FontSize=textboxfontsize*0.6;
textCDf.LineWidth=linewidth;

textCDm=annotation('textbox',[dsx-0.1 dsy+0.08 0.1 0.1],'String',{'$H < 0$','$\Omega_A = \sqrt[\alpha + 1]{\gamma}\Omega$'},'interpreter','latex');
textCDm.EdgeColor=dScolor;
textCDm.FontSize=textboxfontsize;
textCDm.LineWidth=linewidth;

textCDCD=annotation('textbox',[dsx-0.29 dsy-0.06 0.1 0.1],'String',{'$H < 0$','k=+1' '$\Omega_\Lambda + \Omega_A > (\Omega_\Lambda + \Omega_A)_c$'},'interpreter','latex');
textCDCD.EdgeColor=dScolor;
textCDCD.FontSize=textboxfontsize*.7;
textCDCD.LineWidth=linewidth;

textCDe=annotation('textbox',[dsx-0.11 dsy-0.2 0.1 0.1],'String',{'$H < 0$','k=+1' '$\Omega_\Lambda + \Omega_A = (\Omega_\Lambda + \Omega_A)_c$'},'interpreter','latex');
textCDe.EdgeColor=dScolor;
textCDe.FontSize=textboxfontsize*.7;
textCDe.LineWidth=linewidth;


%E and around
txtE=annotation('textbox',[ex-0.05 ey-0.1 0.1 0.1],'String','\textbf{E}','EdgeColor','none','HorizontalAlignment',"center",'interpreter','latex');
txtE.FontSize=fontsize;
txtE.Color=Ecolor;

textEf=annotation('textbox',[ex+0.045 ey-0.02 0.1 0.1],'String','$H < 0$','interpreter','latex');
textEf.EdgeColor=Ecolor;
textEf.FontSize=textboxfontsize;
textEf.LineWidth=linewidth;

textEds=annotation('textbox',[ex-0.12 ey-0.02 0.1 0.1],'String','$H > 0$','interpreter','latex');
textEds.EdgeColor=Ecolor;
textEds.FontSize=textboxfontsize;
textEds.LineWidth=linewidth;



% %arrows

%m arrows
mfx=[mx mx+(fx-mx)/5];
mfy=[my my-(my-fy)/5];
mfarrow=annotation('arrow',mfx,mfy);

mdsx=[mx mx-(mx-dsx)/5];
mdsy=[my my-(my-dsy)/5];
mdsarrow=annotation('arrow',mdsx,mdsy);



%f arrows
fmx=[fx fx-(fx-mx)/5];
fmy=[fy fy+(my-fy)/5];
fmarrow=annotation('arrow',fmx,fmy);

fdsx=[fx fx-0.13];
fdsy=[fy dsy];
fdsarrow=annotation('arrow',fdsx,fdsy);

fex=[fx fx-(fx-ex)/5];
fey=[fy fy-(fy-ey)/5];
fearrow=annotation('arrow',fex,fey);


    %F circular arrow
    space=0.1;
    scale=0.05;
    
    for theta=space:space:3*pi/4
    
    annotation('line',fx+0.03+scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],fy-scale*[sin(theta) sin(theta-space)]);
    annotation('line',fx+0.03+scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],fy+scale*[sin(theta) sin(theta-space)]);
    end
    
    fcar=annotation('arrow',[fx+0.029+scale/sqrt(2)+scale fx+0.029+scale/sqrt(2)+scale],[fy+0.005 fy-0.015]);
    fcar.HeadStyle='plain';




%ds arrows
dsfx=[dsx dsx+0.13];
dsfy=[dsy fy];
dsfarrow=annotation('arrow',dsfx,dsfy);

dsmx=[dsx dsx+(mx-dsx)/5];
dsmy=[dsy dsy+(my-dsy)/5];
dsmarrow=annotation('arrow',dsmx,dsmy);

dsex=[dsx dsx+(ex-dsx)/5];
dsey=[dsy dsy-(dsy-ey)/5];
dsearrow=annotation('arrow',dsex,dsey);


    %dS Circular Arrow
    start=pi/4;
    
    for theta=start:space:pi
    
    annotation('line',dsx-0.04-scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],dsy-scale*[sin(theta) sin(theta-space)]);
    annotation('line',dsx-0.04-scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],dsy+scale*[sin(theta) sin(theta-space)]);
    end
    
    dscar=annotation('arrow',[dsx-0.0385-scale/sqrt(2)-scale dsx-0.0385-scale/sqrt(2)-scale],[dsy+0.005 dsy-0.015]);
    dscar.HeadStyle='plain';


%E arrows

efx=[ex ex+(fx-ex)/5];
efy=[ey ey+(fy-ey)/5];
efarrow=annotation('arrow',efx,efy);


edsx=[ex ex-(ex-dsx)/5];
edsy=[ey ey+(dsy-ey)/5];
edsarrow=annotation('arrow',edsx,edsy);

%maximize the figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])




%     figname1 = strcat('23gamma2chart','.jpg');
% 
%     h=gcf;
%     set(h,'PaperOrientation','landscape');
%     print(figname1,'-djpeg', '-r500');