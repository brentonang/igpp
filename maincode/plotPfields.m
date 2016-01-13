
clear

rc=[0:1/32:1 ones(size(1:32))];    
gc=[0:1/32:1 31/32:-1/32:0];
bc=[ones(size(1:33)) 31/32:-1/32:0];
rgbcolor=[rc;gc;bc]';  
scl=10
color=jet; color(1,:) = 1;

grav = 9.8;
% this scaling is for an isothermal atmosphere (ie only approximate)
aa = -1*grav/(2*287*300);	% factor of 2 for square root
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmaxold = 999;
cleft = 295;            % speed of window at left after Tdelay
fid=fopen('topo.in');
if (fid>0)
    topo = fscanf(fid,'%f',[2,inf])';       % equiv to load topo.in
    iftopo = 1;
    fclose(fid);
end

% ----------------------read and plot pressure perturbations -------

nn = 100;
for jj = 0:nn
         
    clear dpress xvec zvec 
	NT = jj*1000;
	fnamer = ['pr',sprintf('%2.2d',round(NT)),'.bin']
	fid=fopen(fnamer,'r'); if (fid<0) continue; end
	yp=fread(fid,'double'); fclose(fid);
% nxx x mzz are the maximum dimensions of the array
    mzz = yp(1); nxx = yp(2);
	dx = yp(3); dt = yp(5); zmin = yp(6); NTT=yp(4);
% nx x mz are the maximum dimensions of the array that are written to file
	nx = yp(7); mz = yp(8);
    xvec = [0:nx-1]*dx/1000; zvec = (zmin+[0:mz-1]*dx)/1000;
    zmax = zvec(end)*1000;

    matrix = reshape(yp,mz+1,nx+1);
    dpress = matrix(2:end,2:end);  
		 
	[m,n]=size(dpress);	
	scale = repmat(exp(aa*zvec*1000),n,1)';
	dpress= dpress./scale; 	 

    Tdelay = ((mzz+83)*dx)/cleft
    iwindoleft = floor(max(0, ((jj-1)*1000*dt-Tdelay)*cleft/dx));
    xleft = iwindoleft*dx/1000;

    if iwindoleft >5
        dpress(:,1:iwindoleft) = eps;
        cleft = 325;
% probably should do some damping within main code
% maybe put in a PML at left once window moves away from boundary
    end

%	cmax=min(cmaxold,max(max(abs(dpress))));
    cmax=max(max(abs(dpress)));

    if jj==0
        cscale = cmax/2;
    else
        cscale = 20/(NT*dt*330)
    end


% plot the main field
    figure(1),clf,axes('position',[0.22 0.15 0.56 0.75])
    imagesc(xvec,zvec,dpress,cscale*[-1 1]),axis('xy'),colorbar('SouthOutside')
    hold on,plot(xleft*[1 1],[zmin zmax]/1000,'g')
    if iftopo
        hold on
        plot(topo(:,1)/1000,topo(:,2)/1000,'k')
    end
    xlabel('Range (km)'), ylabel('Altitude (km)')
	colormap(rgbcolor)
    title(['Scaled pressure field at ',num2str(dt*NT),' seconds']),grid
	text(xvec(end)*0.1,zmin/1000+(zmax-zmin)*0.0001,['max value = ',num2str(cmax,3)])
        
% .........		              
    str = ['print -djpeg90 plotPfield',num2str(NTT),'3D.jpg']
    eval(str)

% work out what the right and left window limits should be 

	pscale = max(abs(dpress)); pmax=max(pscale);
	pscale=pscale/pmax;	
    figure(2),clf,semilogy(xvec,pscale),ylim([0.001 1.1])
%
	pind =find(pscale>1/100);
	ileft = pind(1);		% because max(pscale)=1
	iright = pind(end);
% approximate window for next 1000 iterations (needs refinement)
	ileft = max(1,ileft-500);
	iright = min(nxx,iright+500);
	hold on,plot(ileft*[1 1]*dx/1000,[0.0001 1],'r')
	hold on,plot(iright*[1 1]*dx/1000,[0.0001 1],'r')
    plot(xleft*[1 1],[0.0001 1],'g')
    text(xleft+1,0.5,'windows based on 1/100 scaling')
    xlabel('Range (km)')
    xlim([0 nx*dx/1000])

	grid
%	[ileft iright]
	str = ['print -djpeg90 plotPmax',num2str(NTT),'.jpg']
 %   if jj>14
 %       eval(str)
 %   end
    cmaxold = cmax;
        
end
