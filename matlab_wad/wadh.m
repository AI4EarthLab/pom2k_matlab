function [art,aru,arv,fsm,dum,dvm,wetmask,h]=wadh(art,aru,arv,fsm,dum,dvm,wetmask, ...
	im,jm,dx,dy,h,hc,hco,hhi,hmin,nwad,slmax,time)
	%{
	% **********************************************************************
	% *                                                                    *
	% * function    :  sets up wad definitions of topography [see fig.1 of *
	% *                o2006]                                              *
	% *                                                                    *
	% **********************************************************************
	%                                                                      !
	%     input (through common block) is "h" defined as follows:          !
	%                                                                      !
	%     "h" is wrt msl (mean sea level) at which level h=0, and h<0 for  !
	%     depth below the msl, and h>0 for "land" above the msl (i.e. for  !
	%     marshes, wetlands, hills & mountains etc.)                       !
	%                                                                      !
	%     also, nwad, hhi, slmax                                           !
	%                                                                      !
	%     outputs:                                                         !
	%                                                                      !
	%     h --> (i.e. becomes) h_old_pom + hhi (which can be zero)         !
	%       i.e. "h" is wrt to some high land-datum (absolute land boundary!
	%       or alb), "hhi" meters above the msl. all cells other than the  !
	%       alb cells are either always wet or can potentially become wet. !
	%                                                                      !
	%     also, fsm, dum, dvm, wetmask, cell areas etc, and hmax           !
	%       see more detailed definitions inside the subroutine            !
	%                                                                      !
	% ---------------------------------------------------------------------!
	%}
	%      integer npos,nneg
	%      real hwatmin
	%      integer i,j
	%     do a rudimentary check on "h:"                                   !                                                                      !
	npos=0; nneg=0;
	for j=1:jm
		for i=1:im
			if(h(i,j)>=0.0)
				npos=+1;
			else
				nneg=-1;
			end
			if(npos*nneg< 0)
				 break;
			end
		end
	end

	if (npos*nneg >= 0)
		fprintf(' stopped in subr. wadh; incorrect defn of h:\n');
		fprintf(' h is one-sign only; see comments in subr.wadh\n');
		fprintf(' npos,nneg =%d %d', npos,nneg);
		% call prxy('undisturbed water depth, h',time,h,im,1,jm,1,0.0)
		% error('')
	end
	hkeep=h ;% keep original for plots etc
	%                                                                      !
	%     hmin & hc (all +ve) are defined as:                              !
	%                                                                      !
	%     h = hmin ---> absolute land area [never flooded], fsm=0.0;       !
	%     h >= hc  ---> water depth wrt high water level, fsm=1.0;         !
	%                                                                      !
	%     note that hc is defined in main --                               !
	%     hmin=0.01; hc=0.05; hco=hc       % hco is used to avoid round-
	hmin=0.01; hco=hc;               % hco is used to avoid round-
	hc=hc*1.0001;                    % off in initial elevation
	%     hc=hc*1.01                       % tne:!wad:- using a larger
	%                 multiplying factor, "1.01" instead of "1.0001",      !
	%                 can help eliminate singular wet spots due to roundoff!
	%                                                                      !
	fprintf(' outputs from subr. wadh ---------------------\n');
	fprintf(' hmin,hc,hco =%e %e %e\n', hmin,hc,hco);
	%                                                                      !
	%     reverse sign of "h", and shift reference if hhi.ne.0, see below..!
	%     (note that hhi is defined in main)                               !
	h=hhi-h;

	if (nwad==1)
		%       if nwad=1, then hhi.ne.0.0, and line "h(:,:)=-h(:,:)+hhi" above  !
		%       already shifts the reference level from msl to alb, the absolute !
		%       land boundary, according to o2006, fig.1:                        !
		for j=1:jm
			for i=1:im
				fsm(i,j)=1.;	dum(i,j)=1.; dvm(i,j)=1.;
				if(h(i,j)<0.0)              % define absolute land areas (never flood):
					%  temporarily set h=-1.0 (some -ve number)
					h(i,j)=-1.0; fsm(i,j)=0.;
					dum(i,j)=0.; dvm(i,j)=0.;
				end;
			end;
		end;
		for j=1:jm-1;
			for i=1:im
				if(fsm(i,j)==0.&&fsm(i,j+1)~=0.)
					dvm(i,j+1)=0.;
				end
			end;
		end
		for j=1:jm;
			for i=1:im-1
				if(fsm(i,j)==0.&&fsm(i+1,j)~=0.)
					dum(i+1,j)=0.;
				end
			end;
		end
		%{                                                                       !
		%       the above yields the following definitions for h:                !
		%                                                                        !
		%       h=-1, absolute land areas [never flooded], fsm=0.0, always       !
		%       h>=0, wet but potentially dry areas,  fsm also=1.0, always       !
		%                                                                        !
		%       note that slpmax touches fsm=1 points only - wet or potentially  !
		%       wad cells; it changes "h" so a wet cell might become dry; so it  !
		%       must be called here before wetmask is defined:                   !
		%                                                                        !
		%       adjust bottom topography so that cell to cell variations         !
		%       in h do not exceed parameter slmax:                              !
		%}
		if(slmax<1.)
			[h]=slpmax(h,im,jm,fsm,slmax)
		end
		%{
		%       we now define a wetmask, assuming that at t=0, the free surface  !
		%       is at the msl (i.e. elevation=0 wrt msl).  if this run does not  !
		%       begin with time=0 (or if other special free-surface distribution !
		%       is desired), then wetmask needs to be redefined or read in e.g.  !
		%       from a restart file, as done later in main if nread.ne.0         !
		%}
		wetmask=fsm;

		for j=1:jm;
			for i=1:im
				if (h(i,j)<hhi)
					wetmask(i,j)=0.0;
				end
			end;
		end
		%                                                                        !
		%       finalize definitions of "h":                                     !
		%                                                                        !
		%       h=hmin, absolute land areas [never flooded], fsm=0.0, always     !
		%       h>=hc,  wet but potentially dry areas,  fsm also=1.0, always     !
		%                                     ...  but wetmask can be 0 or 1     !
		%                                                                        !
		for j=1:jm;
			for i=1:im
				if (h(i,j)<0.0)
					h(i,j)=hmin;
				else
					h(i,j)=h(i,j)+hc;
				end
			end;
		end
		[art,aru,arv,fsm,dum,dvm]=areas_masks(art,aru,arv,fsm,dum,dvm,...
				im,jm,dx,dy,h,nwad);  % calc. cells' areas etc but do not alter fsm...
		
		% if nwad=1
	else     % nwad.eq.0

		for j=1:jm;
			for i=1:im
				if (h(i,j)<=0.0)
					h(i,j)=hmin;
				end
			end;
		end

		%       set minimum water depth to hwatmin:                              !
		%       hwatmin=min water depth, must be >1 to match subr.areas_masks'   !
		%       definition of land, and large enough to avoid grid cells from    !
		%       becoming dry                                                     !
		%                                                                        !
		hwatmin=10.0;
		for j=1:jm;
			for i=1:im
				if (h(i,j)<hmin&&h(i,j)<hwatmin)
					h(i,j)=hwatmin;
				end
			end;
		end

		[art,aru,arv,fsm,dum,dvm]=areas_masks(art,aru,arv,fsm,dum,dvm,...
		im,jm,dx,dy,h);  % calc. cells' areas etc & fsm if nwad=0
		%       adjust bottom topography so that cell to cell variations         !
		%       in h do not exceed parameter slmax:                              !

		if(slmax<1.e0)
			[h]=slpmax(h,im,jm,fsm,slmax);
		end

		wetmask=fsm; % wetmask is same as fsm if nwad=0

	end   % if (nwad.eq.1) then...else...

	%     print to check:                                                  !

	% call prxy(' input h ',time, hkeep(1,1),im,1,jm,1,0.0)
	% call prxy(' h after wadh ',time, h(1,1),im,1,jm,1,0.0)
	% call prxy(' fsm ',time, fsm(1,1),im,1,jm,1,1.0)
	% call prxy(' wetmask ',time, wetmask(1,1),im,1,jm,1,1.0)
	% tps=fsm-wetmask
	% call prxy('fsm-wetmask',time,tps(1,1),im,1,jm,1,1.0)

	%     double-check masks for nwad=0:                                   !

	if (nwad==0)
		for j=1:jm;
			for i=1:im
				if (fsm(i,j)~=wetmask(i,j))
					fprintf(' stopped, nwad =%d\n',nwad);
					fprintf(' fsm.ne.wetmask @ (i,j) =%d %d\n', i,j);
					error('')
				end
			end;
		end
	end
end
