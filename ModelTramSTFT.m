fs=4; 
tf=300; 
t=0:1/fs:tf-1/fs;
f1=0.04 ;
f2=0.4; 
SLOPE=(f2-f1)./t(end);
F=f1+SLOPE*t;

x=(0.3*cos(2*pi*0.1*t))+(1*cos(2*pi*F.*t));

% % duracion=600; %% TAMAÑO DE LA VENTANA DE ANALISIS EN SEGUNDOS
% % ventana=duracion*fs;
% % 
% %   a=x(1:ventana/2);
% %   aa=flip(a);
% %   b=x(end-(ventana/2):end);
% %   c=flip(b);
% %   
% % X=cat(1,aa,x,c);
% %   
% % gridSegundos=360; %% DESPLAZAMIENTO DEL GRID EN SEGUNDOS
% % puntoscorri=(gridSegundos*4); %+((gridSegundos/2)*4)
% % grid=ventana-puntoscorri;
% % nfft=2*ventana; %% NUMERO DE PUNTOS DE LA STFT

 [S,F,T,P]=spectrogram(x,300,296,300,fs,'yaxis');
% % % 
ii=figure(2)
surf(T,F,P,'FaceColor', 'interp',...
     'EdgeColor','none',...
     'FaceLighting','phong');
 axis tight;
 view(0,90);
 ylabel('Frecuencia (Hz)');
 xlabel('Tiempo (s)');
 title('Representación T-F con STFT del modelo IPFM ');
 colorbar
 colormap('jet');
 ylim([0 0.4]);
% caxis([0 max(max(P(:,:)))/20]);

 print(ii,'IPFMmodelSTFT','-dtiff','-r300');