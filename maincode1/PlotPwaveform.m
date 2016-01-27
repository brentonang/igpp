%   1) read in the waveform data
%   2) get receiver location info
%   3) plot scaled pressure amplitudes (pressure * source distance)

clear

redvel = .350;									% reducing velocity
srcalt=2;
path(path,'/Users/cdh/Catherine/Tmatlab.code')

%*******************************************************************
% ----------------------read and plot  wavforms -------
%*******************************************************************

% linear run,  viscosity included

%NT = 112821;               %for TestC
%NT = 107318;               %for TestB
NT = 100001;                %for TestA
%NT = 77194;                 % for TestX
fname = ['prespCyl',sprintf('%2.2d',round(NT)),'.bin']
fid=fopen(fname,'r'); yp=fread(fid,'double'); fclose(fid);

% mzz x nxx are the maximum dimensions 
% (nxx - max # stations, mzz - max # time steps)
mzz = yp(1); nxx = yp(2); 
% nsta x mz are the maximum dimensions of the array that are written to file
% (actual number of stations and time steps)
nsta = yp(7); mz = yp(8);
dx = yp(3); dt = yp(5); zmin = yp(6); NTT=yp(4);
tvec = [0:mz-2]*dt;

matrix = reshape(yp,mz+1,nsta+1);
recwaves = matrix(2:end-1,2:end);  cmax=max(max(abs(recwaves)));
recvrlocs = matrix(1,2:end);
%ilocs = mod(recvrlocs,100000);
%jlocs = (recvrlocs - ilocs)/100000;

ilocs = floor(recvrlocs);
jlocs = round((recvrlocs - ilocs)*100000);
xvec = (ilocs-0.5)*dx/1000;

% source distance in km
srcdist = sqrt((srcalt - dx/2)^2 + (dx*(ilocs-0.5)).^2)/1000;
rrkm = dx*(ilocs-0.5)/1000;

% 2) copy data in time windows about the signal..................
figure(1),clf,axes('position',[0.1 0.1 0.8 0.8])
nsta = 700;             % only plot to receivers at this distance

for ii = 1:3:nsta
    plot(tvec-srcdist(ii)/redvel,rrkm(ii)+20*srcdist(ii)*recwaves(:,ii),'k')
    hold on
end
axis([-1 40 -1 1+ceil(rrkm(nsta))])
title(['Scaled Pressure Waveforms vs Range from source'])
xlabel(['Reduced Time (s) (Reducing Velocity = ', num2str(redvel*1000),' m/s)'])
%ylabel('Pressure*Source Distance (Pa*km)'),grid
ylabel('Range from source (km)'),grid

orient tall
str = ['print -djpeg90 plotPwaves',num2str(NT),'.jpg']
eval(str)

keyboard

twindo = 43;    % form spectra over this time length
twpre=4;
ntpick = round(twindo/dt);
fixwaves = zeros(ntpick+1,nsta);
predwaves = zeros(ntpick+1,nsta);

% 2) copy data in time windows about the signal..................
% figure(1),clf           %  axes('position',[0.1 0.2 0.8 0.4])
for ii = 1:nsta
%    ii
%    plot(tvec-srcdist(ii)/redvel,rrkm(ii)+3*srcdist(ii)*recwaves(:,ii),'r')
%    hold on
    t1 = max(1,round((srcdist(ii)/redvel)/dt-twpre/dt));
    t2 = t1+ntpick;
%    [ii t1 t2 mz]
%    predwaves(:,ii) = srcdist(ii)*recwaves(t1:t2,ii);
%    plot(tvec(t1:t2)-srcdist(ii)/redvel,rrkm(ii)+3*predwaves(:,ii),'k')
    fixwaves(t1:t2,ii) = recwaves(t1:t2,ii);
% pause(0.1)
end
%xlim([-5-twpre twindo])
%title(['Pressure Waveforms'])
%xlabel('Reduced Time (s) (Reducing Velocity = 355 m/s)')
%ylabel('Pressure * (Source Distance) (kPa*m)'),grid

%tspeclen = 17.6;			% form spectra over this time length
%nspeclen = floor(tspeclen/dt);
fmax=10;

% 3) make a matrix of amplitude spectra vs distance ................
for ii = 1:nsta
    trace = fixwaves(:,ii);
    [ampyy,f] = ampspec(trace,1/dt);
    if (ii == 1)
        ifreq = find(f<=fmax); nfreq = length(ifreq)
        specmat1 = zeros(nfreq,nsta);
    end
    [amax,ind] = max(ampyy);
    fmaxv(ii) = f(ind);
    parsevalEn(ii) = sum(ampyy.*ampyy)*length(trace);
    specmat(:,ii) = 20*log10(ampyy(1:nfreq)+eps);
end
freqv = f(1:nfreq);

figure(3),clf, axes('position',[0.1 0.3 0.7 0.4])
imagesc(rrkm(1:nsta),freqv,specmat+130,[-40 40]), axis xy, grid
colorbar('vert')
% h=gca;set(h,'XTickLabel',[])
ylabel('Frequency (Hz)'), xlabel('Range from source (km)')
title('Spectra as a function of range from source')

str = ['print -djpeg90 SpectraPwaves',num2str(NT),'.jpg']
eval(str)


return


keyboard
[ii t1 t2 ntpts]
xlim([-1 4])
title(['Pressure Waveforms'])
xlabel('Reduced Time (s) (Reducing Velocity = 340 m/s)')
ylabel('Pressure * (Source Distance) (kPa*m)'),grid

return
figure(3),clf,subplot(2,1,1)
plot(xvec(3:end),(max(recwaves(:,3:end))-min(recwaves(:,3:end))).*srcdist(3:end))
grid, ylabel('Scaled waveforms amplitudes')
subplot(2,1,2)
plot(xvec,(jlocs-1)*dx+zmin),grid
ylabel('Altitude (m)'),xlabel('Range (km)')





