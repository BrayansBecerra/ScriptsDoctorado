fs=4;
t = 0:1/fs:300;
x = chirp(t,0.003,300,0.4,'quadratic');

[S,F,T,P]=spectrogram(x,400,396,400,fs,'yaxis');
 
gh=figure(1);
surf(T,F,P./max(P(:,:)),'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');
axis tight;
view(0,90);
ylim([0 0.4]);
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (s)');
%legend('Ventana de analisis de: 32 segundos');
title('Representación T-F de una función Chirp cuadrática con STFT');
%,'FontSize',14,'FontWeight','bold');
colormap('jet');
colorbar
%caxis([0 max(P(:,1))]);
print(gh,'SpectrogramNormalizeChirpSignal','-dtiff','-r300');


%x = 0.4*cos(2*pi*0.04*t)+0.3*cos(2*pi*0.2*t);%+(0.4*rand(1,length(t)));
[S,F,P]=cwt(x,'amor',Fs);
s=abs(S);
ii=figure(2)
  surf(t,F,s, 'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');
axis tight;
view(0,90);
ylim([0 0.4]);
xlim([0 300]);
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (s)');
title('Representación T-F de una función Chirp cuadrática con CWT');
colormap('jet');
colorbar
%print(ii,'ScalogramChirpSignal','-dtiff','-r300');
