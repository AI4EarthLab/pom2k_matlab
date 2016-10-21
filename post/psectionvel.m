function psectionvel(vel,east_e,north_e,east_key,north_key,velkey,isub,scalfac,color)
%
% psectionvel:     plots a horizontal velocity section using arrows
%
% Usage: psectionvel(vel,east_e,north_e,east_key,north_key,velkey,isub,scalfac,[color])
%
% where: vel ..... complex 2-D section of horizontal velocity vectors, with 
%                  components in east and north
%        east_e .. horizontal coordinates of centres of model "cell" (metres)
%        north_e . horizontal coordinates of centres of model "cell" (metres)
%        east_key .. east-position of key arrow
%        north_key . north-position of key arrow
%        velkey .. velocity for key arrow (metres/second)
%        isub .... subsample factor
%        scalfac . scale factor for arrows
%        color ... colour for arrows
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
%
if (nargin<8 | nargin>9)
  help psectionvel;
  return
end
%
if (nargin==8)
  color='black';
end
%
[m,n]=size(vel);
%
vel=vel([isub:isub:m],[isub:isub:n]);
east_e=east_e([isub:isub:m],[isub:isub:n]);
north_e=north_e([isub:isub:m],[isub:isub:n]);
%
ig=find(finite(vel));
%
h=arrows(east_e(ig),north_e(ig),vel(ig),scalfac,color);
%
% Add key:
%
h=arrows(east_key,north_key,velkey,scalfac,color);
text(east_key,north_key-velkey*scalfac,strcat(num2str(velkey),'m/s'));
