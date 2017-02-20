function activeSubplot=plot_mancal_callback(hEdit1,hEdit2,activeSubplot,TS,TM,H)

s1 = get(hEdit1,'string');
s2 = get(hEdit2,'string');


p(1)=str2double(s1);
p(2)=str2double(s2);

if isnan(p(1))
    blink(hEdit1)
    return
end

if isnan(p(2))
    blink(hEdit2)
    return
end



pred = repmat(NaN,size(TS));
IO = TS>0;
pred(IO)=sluginj(p,TS(IO));

if activeSubplot==5
    activeSubplot=activeSubplot+1;
end

subplot(2,4,activeSubplot)

plot(TM,H,'.',TS,pred)
set(gca,'ylim',[0,1],'xlim',[TS(1),TS(end)])
xlabel('Time (h)')
ylabel('Head (m)');
title(['S = ',num2str(p(1)),' ; T = ',num2str(p(2))])

if activeSubplot<8
    activeSubplot=activeSubplot+1;
end



function blink(hEdit)

for u=1:4
    set(hEdit, 'BackgroundColor','r')
    drawnow
    pause(0.25)
    set(hEdit, 'BackgroundColor','w')
    pause(0.25)
    drawnow
end

