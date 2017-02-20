function [hEdit1,hEdit2,hEdit3,hEdit4,hButtonGo,hButtonClf,activeSubplot]=disp_ui(hFig)

uicontrol(hFig,'style','text','position',[5,160,70,20],'string','vegFraction')
uicontrol(hFig,'style','text','position',[5,120,70,20],'string','drainageEff')
uicontrol(hFig,'style','text','position',[5,80,70,20],'string','canopyMax')
uicontrol(hFig,'style','text','position',[5,40,70,20],'string','evapFactor')

hEdit1 = uicontrol(hFig,'style','edit','position',[80,160,40,20],...
    'string','0.6','BackgroundColor','w');
hEdit2 = uicontrol(hFig,'style','edit','position',[80,120,40,20],...
    'string','200','BackgroundColor','w');
hEdit3 = uicontrol(hFig,'style','edit','position',[80,80,40,20],...
    'string','2.5','BackgroundColor','w');
hEdit4 = uicontrol(hFig,'style','edit','position',[80,40,40,20],...
    'string','0.83','BackgroundColor','w');



hButtonGo = uicontrol(hFig,'style','pushbutton','position',[70,10,40,20],...
                     'string','Go',...
                     'Callback',['activeSubplot=plot_mancal_callbac',...
                     'k(hEdit1,hEdit2,hEdit3,hEdit4,activeSubplot);']);
hButtonClf = uicontrol(hFig,'style','pushbutton','position',[20,10,40,20],...
                     'string','Reset',...
                     'Callback', ['figure(hFig),clf,[hEdit1,hEdit2,hEdit3,hEdit4,hButt',...
                     'onGo,hButtonClf,activeSubplot]=disp_ui(hFig);']);
                 

activeSubplot = 1;

set(gcf,'toolbar','figure')