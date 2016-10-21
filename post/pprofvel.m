function pprofvel(vel,height,heightlim,i,j,itime,iplot,ileg)
%
% pprofvel:        plots one velocity profile
%
% Usage: pprofvel(vel,height,heightlim,i,j,itime,iplot,ileg)
%
% where: vel ..... velocity profile (complex)
%        height .. vertical coordinate
%        heightlim .. limits of height for plotting
%        i ....... i-index of cell
%        j ....... j-index of cell
%        itime ... the time index (first output = 1, etc.)
%        iplot ... subplot number
%        ileg .... 1 for legend, otherwise zero
%
% Initial version, JRH 11/12/2001
% itime added to title, JRH 02/01/2002
%
hold on;
grid on;
box on;
%
plot(real(vel),height,'rx-');plot(imag(vel),height,'b+-');
%
set(gca,'YLim',heightlim);
plot([0 0],heightlim,'k -');
%
tit=sprintf('Cell %s (%d, %d); time index %d',char(iplot+64),i,j,itime);
title(tit);
xlabel('Velocity (m/s)');
ylabel('Height (m)');
if ileg == 1
  legend('Eastward','Northward');
end
%
hold off;
