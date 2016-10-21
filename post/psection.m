function psection(var,px,py,ar,tit,cax,label)
%
% psection:        plots a section (an extended version of pcolor)
%
% Usage: psection(var,px,py,ar,tit,[cax],[label])
%
% where: var ..... array of scalar data
%        px ...... horizontal plotting coordinate
%        py ...... vertical plotting coordinate
%        ar ...... aspect ratio of plot (< 1 compresses horizontal axis)
%        tit ..... the title (a string)
%        cax ..... range of data to plot (autoscales if not supplied)
%        label ... legend for colourbar
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
% Title added, JRH 02/01/2002
% Correction to test of narg, JRH 11/1/2002
%
if ( nargin < 5 | nargin > 7 )
  help psection;
  return
end
%
ind=find(~isnan(var));
if(nargin < 6),
 cax(1)=min(var(ind));
 cax(2)=max(var(ind));
end
var(ind)=max(var(ind),cax(1));
var(ind)=min(var(ind),cax(2));
%
[m,n]=size(px);
%
% If necessary, pad var out to be same size as px (assumed to be same as py):
%
[mvar,nvar]=size(var);
%
if(m > mvar)
  var=[var;zeros(m-mvar,nvar)];
end
%
if(n > nvar)
  var=[var zeros(m,n-nvar)];
end
%
if(exist('label')==0)
  label='';
end
%
pcolor(px,py,var);colormap('jet');shading('flat');caxis(cax);
%
h=colorbar('horiz');
xltag=get(h,'Xlabel');
set(xltag,'string',label);
%
title(tit);
xlabel('metres');
ylabel('metres');
%
set(h,'Position',[ .28 .1 .44 .05 ]);
%
set ( gca, 'DataAspectRatio', [1 ar 1] );
