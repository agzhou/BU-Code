 handles.gg2 = gg;
handles.sIQHP2 = sIQHP;
handles.FsIQHP2 = FsIQHP;
handles.rNeg2 = rNeg;
handles.rPos2 = rPos;
guidata(hObject, handles);

gg2 = handles.gg2./handles.gg2(1,1);
gg = gg./gg(1,1);

gg2 = handles.gg2;
sIQHP2 = handles.sIQHP2;
FsIQHP2 = handles.FsIQHP2;
rNeg2 = handles.rNeg2;
rPos2 = handles.rPos2;
%% gg
Fig=figure;
set(Fig,'Position',[300 200 400 1000]);
subplot(5,1,2);
plot(handles.tauCoor*1e3,squeeze(real(gg)),'k');hold on;plot(handles.tauCoor*1e3,squeeze(real(gg2)),'-.k');
xlabel('tau [ms]')
grid on
title(['GG-real,', '[', num2str(handles.slt),']'])
subplot(5,1,3);
plot(handles.tauCoor*1e3,squeeze(imag(gg)),'b');hold on;plot(handles.tauCoor*1e3,squeeze(imag(gg2)),'-.b');
xlabel('tau [ms]')
grid on
title(['GG-imag'])
subplot(5,1,5);
plot(squeeze((gg)),'k');hold on;plot(squeeze((gg2)),'-.k');
xlim([-1 1]);ylim([-1 1])
grid on
title(['GG-complex'])
subplot(5,1,1);
plot(handles.tauCoor*1e3,squeeze(abs(gg)),'r');hold on;plot(handles.tauCoor*1e3,squeeze(abs(gg2)),'-.r');
xlabel('tau [ms]')
ylim([0 1])
grid on
title(['GG-magnitude'])
legend({'stim','base'})
subplot(5,1,4);
plot(handles.tauCoor*1e3,squeeze(angle(gg)),'k');hold on;plot(handles.tauCoor*1e3,squeeze(angle(gg2)),'-.k');
xlabel('tau [ms]'); ylim([-pi, pi]);
grid on
title(['GG-phase'])

%% Fgg
fCoor = linspace(-2500,2500,100);
Fgg = fftshift(fft(gg));
Fgg2 = fftshift(fft(gg2));
Fig=figure;
set(Fig,'Position',[300 200 400 1000]);
subplot(5,1,2);
plot(fCoor,squeeze(real(Fgg)),'k');hold on;plot(fCoor,squeeze(real(Fgg2)),'-.k');
xlabel('tau [ms]')
grid on
title(['GG-real,', '[', num2str(handles.slt),']'])
subplot(5,1,3);
plot(fCoor,squeeze(imag(Fgg)),'b');hold on;plot(fCoor,squeeze(imag(Fgg2)),'-.b');
xlabel('tau [ms]')
grid on
title(['GG-imag'])
subplot(5,1,5);
plot(squeeze((Fgg)),'k');hold on;plot(squeeze((Fgg2)),'-.k');
xlim([-1 1]);ylim([-1 1])
grid on
title(['GG-complex'])
subplot(5,1,1);
plot(fCoor,squeeze(abs(Fgg)),'r');hold on;plot(fCoor,squeeze(abs(Fgg2)),'-.r');
xlabel('tau [ms]')
grid on
title(['GG-magnitude'])
legend({'stim','base'})
subplot(5,1,4);
plot(fCoor,squeeze(angle(Fgg)),'k');hold on;plot(fCoor,squeeze(angle(Fgg2)),'-.k');
xlabel('tau [ms]'); ylim([-pi, pi]);
grid on
title(['GG-phase'])

%% fft(gg) for abs real imag angle
fCoor = linspace(-2500,2500,100);
Fig=figure;
set(Fig,'Position',[300 200 400 1000]);
subplot(5,1,2);
plot(fCoor,abs(fftshift(fft((real(gg))))),'k');hold on;plot(fCoor,abs(fftshift(fft((real(gg2))))),'-.k');
xlabel('tau [ms]')
grid on
title(['GG-real,', '[', num2str(handles.slt),']'])
subplot(5,1,3);
plot(fCoor,abs(fftshift(fft((imag(gg))))),'b');hold on;plot(fCoor,abs(fftshift(fft((imag(gg2))))),'-.b');
xlabel('tau [ms]')
grid on
title(['GG-imag'])
subplot(5,1,5);
plot(squeeze((Fgg)),'k');hold on;plot(squeeze((Fgg2)),'-.k');
xlim([-1 1]);ylim([-1 1])
grid on
title(['GG-complex'])
subplot(5,1,1);
plot(fCoor,abs(fftshift(fft((abs(gg))))),'r');hold on;plot(fCoor,abs(fftshift(fft((abs(gg2))))),'-.r');
xlabel('tau [ms]')
grid on
title(['GG-magnitude'])
legend({'stim','base'})
subplot(5,1,4);
plot(fCoor,abs(fftshift(fft((angle(gg))))),'k');hold on;plot(fCoor,abs(fftshift(fft((angle(gg2))))),'-.k');
xlabel('tau [ms]'); 
grid on
title(['GG-phase'])

%% sIQ and FsIQ
fig=figure;
set(fig,'Position',[200 400 1000 700])
subplot(211);plot(tCoor,squeeze(sIQHP),'r');hold on; plot(tCoor, squeeze(sIQHP2),'b'); 
title('sIQHP'); 
xlabel('time lag, [ms]')
legend({'stim','base'})

fCoor = linspace(-2500,2500,1000);
% fig=figure;
% set(fig,'Position',[200 400 1000 350])
subplot(212);plot(fCoor,squeeze(abs(FsIQHP)),'r');hold on; plot(fCoor, squeeze(abs(FsIQHP2)),'b'); 
title('FsIQHP'); 
xlabel('Frequency, [Hz]')
legend({['stim: Pneg=', num2str(rNeg),', Ppos=',num2str(rPos)],['base: Pneg=', num2str(rNeg2),', Ppos=',num2str(rPos2)]});