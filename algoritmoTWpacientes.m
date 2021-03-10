%%load 'RRrem_p1s1.txt';
%'TacRRremSpli4Hz.txt'

fs=4;
%t = 0:1/fs:600-(1/fs);
%x = chirp(t,0.4,599,0.003,'quadratic');


nuevo=load ('RRrem_ST1s1.txt');  %%corregir linea de tiempo del archivo 12
Horas=input('Ingrese el tiempo de terapia(horas): ')
Minutos=input('Ingrese el delta de tiempo (minutos): ')
Soloterapia=nuevo(Minutos*60*4:(Minutos*60*4)+(Horas*60)*60*4,:);

%tamano=length(nuevo(:,1));
t=(Soloterapia(:,1)-(27*60));
x=Soloterapia(:,2); %(Aqui ya estan dado en segundos)

figure(1)
plot(nuevo(:,1)/60,nuevo(:,2));
hold on
plot(Soloterapia(:,1)/60,Soloterapia(:,2));
%hold on
%plot((Soloterapia(:,1)-(27*60))/60,Soloterapia(:,2));

%length(x);
y=x;
%y=detrend(x);
T=t';

figure(2)
plot(t,x,'r')
xlabel('tiempo (s)')
ylabel('RR')
title('Tacograma real')
%figure
%cwt(y,'morse',fs);
hh=figure(3)
cwt(y,'amor',fs);
colormap('jet');
colorbar
caxis([0 0.042]);
%print(gcf,'-dtiff',strcat('C:\imagenesTW','ScalogramaP1.tif'));
%print(hh,'ScalogramaP11','-dtiff','-r300');

%saveas(gcf,'ScalogramaP1.tif');

[S,F,P]=cwt(y,'amor',fs);
%[wt, f]=cwt(y,'morse',fs);
%WT=real(wt);
s=abs(S);
ii=figure(4)
% % 
  surf(T/60,F,s, 'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');
           
axis tight;
% % % axis ij;
% % % view(-45,35);
view(0,90);
ylim([0 0.4]);
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (min)');
title('Representacion tiempo-frecuencia de la HRV con CWT');
% % % %zlim([0 0.08]);
% % % title('Representacion tiempo-frecuencia con CWT');
colormap('jet');
colorbar
% % % caxis([0 max(max(WT(:,:)))]);
caxis([0 0.042]);
print(ii,'ScalogramaP12','-dtiff','-r300');

[cfsVLF,f] = cwt(y,fs);
[cfsLF,f] = cwt(y,fs);
[cfsHF,f] = cwt(y,fs);
%[cfs2,f] = cwt(y,fs);
%T1 = 0;  %T2 = t(end,1);
F1 = 0.003;   F2 = 0.04;
F3=0.15; F4=0.4; 
%cfs(f < F1 & f > F2, t> T1 & t < T2) = 0;


cfsVLF(f < F1, t>0) = 0;
cfsVLF(f>F2, t>0)=0;

cfsLF(f < F2, t>0) = 0;
cfsLF(f>F3, t>0)=0;

cfsHF(f < F3, t>0) = 0;
cfsHF(f>F4, t>0)=0;

xrecVLF = icwt(cfsVLF);
xrecLF = icwt(cfsLF);
xrecHF = icwt(cfsHF);

jj=figure(5);
subplot(3,2,1);
plot(t/60,xrecVLF);
xlabel('Tiempo (min)')
ylabel('RR (s)')
title('Señal reconstruida por ICWT con componentes en la banda de VLF');
xlim([0 180]);

subplot(3,2,2);
[S,F,P]=cwt(xrecVLF,fs);
s=abs(S);
surf(T/60,F,s, 'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');     
axis tight;
view(0,90);
ylim([0 0.4]);
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (min)');
title('Representacion tiempo-frecuencia de la banda de VLF con CWT');
colormap('jet');
colorbar
caxis([0 0.042]);

subplot(3,2,3);
plot(t/60,xrecLF);
xlabel('Tiempo (min)')
ylabel('RR (s)')
title('Señal reconstruida por ICWT con componentes en la banda de LF');
xlim([0 180]);

subplot(3,2,4);
[S,F,P]=cwt(xrecLF,fs);
s=abs(S);
surf(T/60,F,s, 'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');     
axis tight;
view(0,90);
ylim([0 0.4]);
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (min)');
title('Representacion tiempo-frecuencia de la banda de LF con CWT');
colormap('jet');
colorbar
caxis([0 0.042]);

subplot(3,2,5);
plot(t/60,xrecHF);
xlabel('Tiempo (min)')
ylabel('RR (s)')
title('Señal reconstruida por ICWT con componentes en la banda de HF');
xlim([0 180]);

subplot(3,2,6);
[S,F,P]=cwt(xrecHF,fs);
s=abs(S);
surf(T/60,F,s, 'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');     
axis tight;
view(0,90);
ylim([0 0.4]);
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (min)');
title('Representacion tiempo-frecuencia de la banda de HF con CWT');
colormap('jet');
colorbar
caxis([0 0.042]);
print(jj,'icwtScalogramaP1','-dtiff','-r300');


figure
cwt(xrecLF,fs);
colormap('jet');

figure
cwt(xrecHF,fs);
colormap('jet');

figure
subplot(5,1,1);
plot(t,x,'r');
grid on;
title('Original Signal');
subplot(5,1,2);
plot(t,detrend(y));
hold on
plot(t,(xrecVLF+xrecLF+xrecHF),'red')
grid on;
title('Original Signal WITHOUT TREND');

subplot(5,1,3);
plot(t,xrecVLF)
hold on
grid on;
title('Signal in VLF');

%bandpower(xrecVLF)

% subplot(5,1,4);
% plot(t,xrecLF)
% hold on
% grid on;
% title('Signal in LF');
% 
% subplot(5,1,5);
% plot(t,xrecHF)
% hold on
% grid on;
% title('Signal in HF');

