Fs=4; 
tf=300; 
t=0:1/Fs:tf-1/Fs;
f1=0.04 ;
f2=0.4; 
SLOPE=(f2-f1)./t(end);
F=f1+SLOPE*t;


x=(0.3*cos(2*pi*0.1*t))+(1*cos(2*pi*F.*t));

jj=figure(1)
plot(t,x)
xlabel('tiempo (s)');
ylabel('Amplitud (u.a.)')
title('Modelo IPFM x=(0.3*cos(2*pi*0.1*t))+(1*cos(2*pi*Fv.*t))','FontSize',14,'FontWeight','bold');
print(jj,'ScalogramaModelSignal','-dtiff','-r300');

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
title('Representacion T-F con CWT del modelo IPFM');
colormap('jet');
colorbar
print(ii,'ScalogramaModelSignal','-dtiff','-r300');
%colorbar
%caxis([0 0.042]);

