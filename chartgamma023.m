close



Fcolor=[0.5 0 1];
dScolor=.6*[0 0.7 0.5];
Mcolor=[0.7 0.2 0.2];


fontsize=20;
linewidth=1.5;

mx=0.5;
my=0.9;
fx=0.7;
fy=0.3;
dsx=1-fx;
dsy=fy;


%lines
annotation('line',[mx fx],[my fy]);
annotation('line',[fx dsx],[fy dsy]);
annotation('line',[dsx mx],[dsy my]);


%texts


%M and around
txtM=annotation('textbox',[mx-0.05 my-0.02 0.1 0.1],'String','\textbf{M}','EdgeColor','none','HorizontalAlignment',"center",'interpreter','latex');
txtM.FontSize=fontsize;
txtM.Color=Mcolor;

textMf=annotation('textbox',[mx+0.06 my-0.18 0.1 0.1],'String','$H > 0, \Omega_\Lambda = 0, \Omega_A=0$','interpreter','latex');
textMf.EdgeColor=Mcolor;
textMf.LineWidth=linewidth;

textMds=annotation('textbox',[mx-0.31 my-0.18 0.1 0.1],'String',{'$H > 0$', '$\Omega_\Lambda \neq 0$ or $\Omega_A \neq 0$'},'interpreter','latex');
textMds.EdgeColor=Mcolor;
textMds.LineWidth=linewidth;



%F and around
txtF=annotation('textbox',[fx+0.01 fy-0.06 0.1 0.1],'String','\textbf{F}','EdgeColor','none','interpreter','latex');
txtF.FontSize=.7*fontsize;
txtF.Color=Fcolor;

textFm=annotation('textbox',[fx-0.01 fy+0.11 0.1 0.1],'String',{'$H < 0$, k = -1', '$\Omega_A=0, \Omega_\Lambda=0$'},'interpreter','latex');
textFm.EdgeColor=Fcolor;
textFm.LineWidth=linewidth;

textFds=annotation('textbox',[fx-0.18 fy-0.14 0.1 0.1],'String',{'$H > 0$, k = 0'},'interpreter','latex');
textFds.EdgeColor=Fcolor;
textFds.LineWidth=linewidth;

textFf=annotation('textbox',[fx+0.13 fy-0.06 0.1 0.1],'String',{'$H < 0$', 'k = +1'},'interpreter','latex');
textFf.EdgeColor=Fcolor;
textFf.LineWidth=linewidth;





txtdS=annotation('textbox',[dsx-0.08 dsy-0.06 0.1 0.1],'String','\textbf{CD}','EdgeColor','none','interpreter','latex');
txtdS.FontSize=.7*fontsize;
txtdS.Color=dScolor;

textDsf=annotation('textbox',[dsx+0.02 dsy-0.14 0.1 0.1],'String',{'$H < 0$, k = 0'},'interpreter','latex');
textDsf.EdgeColor=dScolor;
textDsf.LineWidth=linewidth;

textDsm=annotation('textbox',[dsx-0.16 dsy+0.08 0.1 0.1],'String','$H < 0$, k = -1','interpreter','latex');
textDsm.EdgeColor=dScolor;
textDsm.LineWidth=linewidth;

textDsds=annotation('textbox',[dsx-0.25 dsy-0.07 0.1 0.1],'String',{'$H < 0$', 'k = +1'},'interpreter','latex');
textDsds.EdgeColor=dScolor;
textDsds.LineWidth=linewidth;



% %arrows


mfx=[mx mx+(fx-mx)/5];
mfy=[my my-(my-fy)/5];
mfarrow=annotation('arrow',mfx,mfy);


fmx=[fx fx-(fx-mx)/5];
fmy=[fy fy+(my-fy)/5];
fmarrow=annotation('arrow',fmx,fmy);


fdsx=[fx fx-0.07];
fdsy=[fy dsy];
fdsarrow=annotation('arrow',fdsx,fdsy);


dsfx=[dsx dsx+0.07];
dsfy=[dsy fy];
dsfarrow=annotation('arrow',dsfx,dsfy);



dsmx=[dsx dsx+(mx-dsx)/5];
dsmy=[dsy dsy+(my-dsy)/5];
dsmarrow=annotation('arrow',dsmx,dsmy);



mdsx=[mx mx-(mx-dsx)/5];
mdsy=[my my-(my-dsy)/5];
mdsarrow=annotation('arrow',mdsx,mdsy);

%F circular arrow
space=0.1;
scale=0.05;

for theta=space:space:3*pi/4

    annotation('line',fx+0.03+scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],fy-scale*[sin(theta) sin(theta-space)]);
    annotation('line',fx+0.03+scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],fy+scale*[sin(theta) sin(theta-space)]);
end

fcar=annotation('arrow',[fx+0.03+scale/sqrt(2)+scale fx+0.03+scale/sqrt(2)+scale],[0.295 0.285]);
fcar.HeadStyle='plain';


%dS Circular Arrow
start=pi/4;

for theta=start:space:pi

    annotation('line',dsx-0.04-scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],dsy-scale*[sin(theta) sin(theta-space)]);
    annotation('line',dsx-0.04-scale/sqrt(2)+scale*[cos(theta) cos(theta-space)],dsy+scale*[sin(theta) sin(theta-space)]);
end

dscar=annotation('arrow',[dsx-0.04-scale/sqrt(2)-scale dsx-0.04-scale/sqrt(2)-scale],[0.295 0.285]);
dscar.HeadStyle='plain';

    
    figname1 = strcat('gamma023chart','.jpg');

    h=gcf;
    set(h,'PaperOrientation','landscape');
    print(figname1,'-djpeg', '-r500');