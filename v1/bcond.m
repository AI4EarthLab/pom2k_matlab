function [elf,uaf,vaf,uf,vf,w] = bcond(idx,elf,uaf,vaf,uf,vf,w,...
                                 im,jm,kb,imm1,jmm1,kbm1,...
                                 fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                                 dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Applies open boundary conditions.                   *
% *                                                                    *
% *                Closed boundary conditions are automatically        *
% *                enabled through specification of the masks, dum,    *
% *                dvm and fsm, in which case the open boundary        *
% *                conditions, included below, will be overwritten.    *
% *                                                                    *
% *                                The C-Grid                          *
% *                                **********                          *
% *                                                                    *
% *                The diagram is for the case where u and v are the   *
% *                primary boundary conditions together with t and     *
% *                s (co-located with el)                              *
% *                                                                    *
% *                All interpolations are centered in space except     *
% *                those at lateral open boundary where an upstream    *
% *                Horizontal locations of e(el), t and s (etc.) are   *
% *                coincident.                                         *
% *                                                                    *
% *                People not acquainted with sigma coordinates have   *
% *                often asked what kind of boundary condition is      *
% *                applied along closed horizontal boundaries.         *
% *                Although the issue is not as important as it might  *
% *                be  for z-level grids, a direct answer is "half-    *
% *                slip" which, of course, is between free slip and    *
% *                non-slip.                                           *
%
%
% East and West end points for the C-grid in POM.
%
%                      west
%
%           v(1,j+1)=0           v(2,j+1) . . . .
%
%     ----<---<----<-----
%     |                 |
%  u(1,j)   el(1,j)   u(2,j)=BC  el(2,j)   u(3,j) . . . .
%             |                   |
%             -----<----<----<-----
%
%           v(1,j)=0              v(2,j) . . . .
%
%                                                    east
%
%                              . . . .  v(im-1,j+1)           v(im,j+1)=0
%
%
%                 . . .  .  u(im-1,j)   el(im-1,j)  u(im,j)=BC  el(im,j)
%                                            |                   |
%                                            ----->----->---->----
%
%                              . . . .   v(im-1,j)             v(im,j)=0
%
%  Notes:
%    1. The suffixes, f  or af, have been deleted.
%    2. All variables NOT designated as boundary condition (=BC) or set to
% zero or obtained from an interior point are calculated points.
%    3. u(1,j) is never used but is obtained from the interior point for
% cosmetic output. Its counterpart, u(im+1,j), fores not exist.
%    4. v=0 at i=1 and i=im are used as open inflow BC s unless specified
% otherwise.
%    5. The south and north extremal points are obtained from the above by
% permuting u to v, v to u, i to j and j to i.


% **********************************************************************
global fsm_3d

if(idx==1)
    %
    %-----------------------------------------------------------------------
    %
    %     External (2-D) boundary conditions:
    %
    %     In this example, the governing boundary conditions are a radiation
    %     condition on uaf in the east and in the west, and vaf in the north
    %     and south. The tangential velocities are set to zero on both
    %     boundaries. These are only one set of possibilities and may not
    %     represent a choice which yields the most physically realistic
    %     result.
    %
    %     Elevation (in this application, elevation is not a primary
    %     boundary condition):
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %     |0 0 ... 0 0|      |0 1 ... 0 0|
    %     |1 1 ... 0 0|      |0 1 ... 0 0|
    % R = |0 0 1 ... 0|  L = |0 0 1 ... 0|
    %     |0 0 ... 1 1|      |0 0 ... 1 0|
    %     |0 0 ... 0 0|      |0 0 ... 1 0|
    %%%%%%%%%%%%%%%%%%%%%%%%
%     L = eye(im);
%     R = eye(jm);
%     elf = (L * elf) * R;
%     elf = elf .* fsm;   
    elf(1, :) = elf(2, :);
    elf(im, :) = elf(imm1, :);
    elf(:, 1) = elf(:, 2);
    elf(:, jm) = elf(:, jmm1);
    elf = elf .* fsm; 
    return
    
elseif(idx==2)
    %
    %     External (2-D) velocity:
      tmph = h;
      tmpel = el;
      tmph(im, :) = tmph(imm1, :);
      tmpel(im, :) = tmpel(imm1, :);
      uab1 = repmat(uabe, im, 1);
      el1 = repmat(ele, im, 1);
      tmpuaf = uab1 + rfe .* sqrt(grav ./(tmph)) .* (tmpel - el1);
      tmpuaf = ramp .* tmpuaf;
      uaf(im, 2:jmm1) = tmpuaf(im, 2:jmm1);
     
      
      tmph(1, :) = tmph(2, :);
      tmpel(1, :) = tmpel(2, :);
      uab2 = repmat(uabw, im, 1);
      el2 = repmat(elw, im, 1);
      tmpuaf = uab2 - rfw .* sqrt(grav ./(tmph)) .* (tmpel - el2);
      tmpuaf = ramp .* tmpuaf;
      uaf(1, 2:jmm1) = tmpuaf(1, 2:jmm1);
      uaf(2, 2:jmm1) = tmpuaf(1, 2:jmm1);
      
      vaf(im, 2:jmm1) = 0.e0;
      vaf(1, 2:jmm1) = 0.e0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      tmph = h;
      tmpel = el;
      tmph(:, jm) = tmph(:, jmm1);
      tmpel(:, jm) = tmpel(:, jmm1);
      vab1 = repmat(vabn, jm, 1);
      vab1 = vab1';
      el3 = repmat(eln, jm, 1);
      el3 = el3';
      tmpvaf = vab1 + rfn .* sqrt(grav./(tmph)).*(tmpel - el3);
      tmpvaf = ramp .* tmpvaf;
      vaf(2:imm1, jm) = tmpvaf(2:imm1, jm);
      
      tmph(:, 1) = tmph(:, 2);
      tmpel(:, 1) = tmpel(:, 2);
      vab2 = repmat(vabs, jm, 1);
      vab2 = vab2';
      el4 = repmat(els, jm, 1);
      el4 = el4';
      tmpvaf = vab2 - rfs .* sqrt(grav./(tmph)).*(tmpel - el4);
      tmpvaf = ramp .* tmpvaf;
      vaf(2:imm1, 1) = tmpvaf(2:imm1, 1);
      vaf(2:imm1, 2) = tmpvaf(2:imm1, 1);
      
      uaf(2:imm1, jm) = 0.e0;
      uaf(2:imm1, 1) = 0.e0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uaf(:,:) = uaf(:,:) .* dum(:,:);
    vaf(:,:) = vaf(:,:) .* dvm(:,:);
    return
    %
elseif(idx==3)
    %
    %-----------------------------------------------------------------------
    %
    %     Internal (3-D) boundary conditions:
    %
    %     Velocity (radiation conditions; smoothing is used in the direction
    %     tangential to the boundaries):
      h_3d = repmat(h, 1, 1, kb);
      tmpu1 = u(imm1, :, :);
      tmpu2 = u(im, :, :);
      tmpga = sqrt(h_3d(im, :, :) ./ hmax);
      
      tmpuf = tmpga .* (AXF1(AXB1(tmpu1))) + (1.e0 - tmpga) .* (AXF1(AXB1(tmpu2)));
      uf(im, 2:jmm1, 1:kbm1) = tmpuf(1, 2:jmm1, 1:kbm1);
      vf(im, 2:jmm1, 1:kbm1) = 0.e0;
      
      tmpu1 = u(3, :, :);
      tmpu2 = u(2, :, :);
      tmpga = sqrt(h_3d(1, :, :) ./ hmax);
      
      tmpuf = tmpga .* (AXF1(AXB1(tmpu1))) + (1.e0 - tmpga) .* (AXF1(AXB1(tmpu2)));
      uf(2, 2:jmm1, 1:kbm1) = tmpuf(1, 2:jmm1, 1:kbm1);
      uf(1, 2:jmm1, 1:kbm1) = tmpuf(1, 2:jmm1, 1:kbm1);
      vf(1, 2:jmm1, 1:kbm1) = 0.e0;
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

      tmpv1 = v(:, jmm1, :);
      tmpv2 = v(:, jm, :);
      tmpga = sqrt(h_3d(:, jm, :) ./ hmax);
      
      tmpvf = tmpga .* (AXF2(AXB2(tmpv1))) + (1.e0 - tmpga) .* (AXF2(AXB2(tmpv2)));
      vf(2:imm1, jm, 1:kbm1) = tmpvf(2:imm1,1,1:kbm1);
      uf(2:imm1, jm, 1:kbm1) = 0.e0;
      
      tmpv1 = v(:, 3, :);
      tmpv2 = v(:, 2, :);
      tmpga = sqrt(h_3d(:, 1, :) ./ hmax);
      
      tmpvf = tmpga .* (AXF2(AXB2(tmpv1))) + (1.e0 - tmpga) .* (AXF2(AXB2(tmpv2)));
      vf(2:imm1, 2, 1:kbm1) = tmpvf(2:imm1,1,1:kbm1);
      vf(2:imm1, 1, 1:kbm1) = tmpvf(2:imm1,1,1:kbm1);
      uf(2:imm1, 1, 1:kbm1) = 0.e0;
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     tmpdum=repmat(dum,1,1,kb);
     tmpdvm=repmat(dvm,1,1,kb);
     
     tmpdum(:, :, kb) = 1.e0;
     tmpdvm(:, :, kb) = 1.e0;
     
     uf = uf .* tmpdum;
     vf = vf .* tmpdvm;

    return
    %
elseif(idx==4)
    %
    %     Temperature and salinity boundary conditions (using uf and vf:
    %     respectively):
    %      
     dx_3d = repmat(dx, 1, 1, kb);
     dt_3d = repmat(dt, 1, 1, kb);
     tmpdtx = dt_3d(imm1, :, :);
     tmpdx1 = dx_3d(im, :, :);
     tmpdx2 = dx_3d(imm1, :, :);  
     tmpzzx = zeros(1, jm, kb);
     tmpzzx(1,:,:) = repmat(zz, jm, 1);
     
     B1x = zeros(1, jm, kb);
     B2x = zeros(1, jm ,kb);
     
     tmpu1x = (2.e0 .* u(im, :, :) .* dti) ./ (tmpdx1 + tmpdx2);
     A1x = tmpu1x;
     A2x = tmpu1x;
     
     A1x(tmpu1x > 0.e0) = 0.e0;
     A1x(tmpu1x <= 0.e0) = 1.e0;
     A2x(tmpu1x > 0.e0) = 1.e0;
     A2x(tmpu1x <= 0.e0) = 0.e0;
     B1x(:,:,2:kb-2) = 1.e0;
     B2x(B1x < .9e0) = 1.e0;

     tmptbe(1,:,:) = tbe;
     tmpsbe(1,:,:) = sbe;
     
     tmpt1x = t(im, :, :) - tmpu1x .* (tmptbe - t(im, :, :));
     tmps1x = s(im, :, :) - tmpu1x .* (tmpsbe - s(im, :, :));
     
     tmpt2x = t(im, :, :) - tmpu1x .* (t(im, :, :) - t(imm1, :, :));
     tmps2x = s(im, :, :) - tmpu1x .* (s(im, :, :) - s(imm1, :, :));
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %       |1  1  0 ...   0  0|
     %       |0  0  1 ...   0  0|
     % R1 =  |0 -1  0 ...   1  0|
     %       |0  0 -1 ...   0  0|
     %       |0  0  0 ...  -1  1|
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     tmpwmx = AZF1(w(imm1, :, :)) .* dti ./ (Rk1(tmpzzx) .* tmpdtx);
     tmpt3x = tmpt2x - tmpwmx .*  (Rk1(t(imm1, :, :)));
     tmps3x = tmps2x - tmpwmx .* (Rk1(s(imm1, :, :)));
     
     tmptx = tmpt1x .* A1x + (tmpt2x .* B2x + B1x .* tmpt3x).*A2x;
     tmpsx = tmps1x .* A1x + (tmps2x .* B2x + B1x .* tmps3x).*A2x;
     
     uf(im, 1:jm, 1:kbm1) = tmptx(1, :, 1:kbm1);
     vf(im, 1:jm, 1:kbm1) = tmpsx(1, :, 1:kbm1);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     tmpdtx = dt_3d(2, :, :);   
     tmpdx1 = dx_3d(1, :, :);
     tmpdx2 = dx_3d(2, :, :);
     
     tmpu1x = (2.e0 .* u(2, :, :) .* dti) ./ (tmpdx1 + tmpdx2);
     A1x = tmpu1x;
     A2x = tmpu1x;
     
     A1x(tmpu1x >= 0.e0) = 1.e0;
     A1x(tmpu1x < 0.e0) = 0.e0;
     A2x(tmpu1x >= 0.e0) = 0.e0;
     A2x(tmpu1x < 0.e0) = 1.e0;
     
     tmptbw(1,:,:) = tbw;
     tmpsbw(1,:,:) = sbw;
     
     tmpt1x = t(1, :, :) - tmpu1x .* (t(1, :, :) - tmptbw);
     tmps1x = s(1, :, :) - tmpu1x .* (s(1, :, :) - tmpsbw);
     
     tmpt2x = t(1, :, :) - tmpu1x .* (t(2, :, :) - t(1, :, :));
     tmps2x = s(1, :, :) - tmpu1x .* (s(2, :, :) - s(1, :, :));
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %       |1  1  0 ...   0  0|
     %       |0  0  1 ...   0  0|
     % R1 =  |0 -1  0 ...   1  0|
     %       |0  0 -1 ...   0  0|
     %       |0  0  0 ...  -1  1|
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     tmpwmx = AZF1(w(2, :, :)) .* dti ./ (Rk1(tmpzzx) .* tmpdtx);
     tmpt3x = tmpt2x - tmpwmx .* (Rk1(t(2, :, :)));
     tmps3x = tmps2x - tmpwmx .* (Rk1(s(2, :, :)));
     
     tmptx = tmpt1x .* A1x + (tmpt2x .* B2x + B1x .* tmpt3x).*A2x;
     tmpsx = tmps1x .* A1x + (tmps2x .* B2x + B1x .* tmps3x).*A2x;
     
     uf(1, 1:jm, 1:kbm1) = tmptx(1, :, 1:kbm1);
     vf(1, 1:jm, 1:kbm1) = tmpsx(1, :, 1:kbm1);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     dy_3d = repmat(dy, 1, 1, kb);
     tmpdty = dt_3d(:, jmm1, :); 
     tmpdy1 = dy_3d(:, jm, :);
     tmpdy2 = dy_3d(:, jmm1, :);
     tmpzzy(:,1,:) = repmat(zz, im, 1);

     B1y = zeros(im, 1, kb);
     B2y = zeros(im, 1, kb);
   
     tmpu1y = (2.e0 .* v(:, jm, :) .* dti) ./ (tmpdy1 + tmpdy2);
     A1y = tmpu1y;
     A2y = tmpu1y;
     
     A1y(tmpu1y > 0.e0) = 0.e0;
     A1y(tmpu1y <= 0.e0) = 1.e0;
     A2y(tmpu1y > 0.e0) = 1.e0;
     A2y(tmpu1y <= 0.e0) = 0.e0;
     B1y(:,:, 2:kb-2) = 1;
     B2y(B1y < .9e0) = 1;

     tmptbn(:,1,:) = tbn;
     tmpsbn(:,1,:) = sbn;
     
     tmpt1y = t(:, jm, :) - tmpu1y .* (tmptbn - t(:, jm, :));
     tmps1y = s(:, jm, :) - tmpu1y .* (tmpsbn - s(:, jm, :));
     
     tmpt2y = t(:, jm, :) - tmpu1y .* (t(:, jm, :) - t(:, jmm1, :));
     tmps2y = s(:, jm, :) - tmpu1y .* (s(:, jm, :) - s(:, jmm1, :));
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %       |1  1  0 ...   0  0|
     %       |0  0  1 ...   0  0|
     % R1 =  |0 -1  0 ...   1  0|
     %       |0  0 -1 ...   0  0|
     %       |0  0  0 ...  -1  1|
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     tmpwmy = AZF(w(:, jmm1, :)) .* dti ./ (Rk2(tmpzzy) .* tmpdty);
     tmpt3y = tmpt2y - tmpwmy .* (Rk2(t(:, jmm1, :)));
     tmps3y = tmps2y - tmpwmy .* (Rk2(s(:, jmm1, :)));
     
     tmpty = tmpt1y .* A1y + (tmpt2y .* B2y + B1y .* tmpt3y).*A2y;
     tmpsy = tmps1y .* A1y + (tmps2y .* B2y + B1y .* tmps3y).*A2y;

     uf(1:im, jm, 1:kbm1) = tmpty(:,1, 1:kbm1);
     vf(1:im, jm, 1:kbm1) = tmpsy(:,1, 1:kbm1);       
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     tmpdty = dt_3d(:, 2, :);     
     tmpdy1 = dy_3d(:, 1, :);
     tmpdy2 = dy_3d(:, 2, :);
     
     tmpu1y = (2.e0 .* v(:, 2, :) .* dti) ./ (tmpdy1 + tmpdy2);
     A1y = tmpu1y;
     A2y = tmpu1y;
     
     A1y(tmpu1y >= 0.e0) = 1.e0;
     A1y(tmpu1y < 0.e0) = 0.e0;
     A2y(tmpu1y >= 0.e0) = 0.e0;
     A2y(tmpu1y < 0.e0) = 1.e0;
     
     tmptbs(:,1,:) = tbs;
     tmpsbs(:,1,:) = sbs;
     
     tmpt1y = t(:, 1, :) - tmpu1y .* (t(:, 1, :) - tmptbs);
     tmps1y = s(:, 1, :) - tmpu1y .* (s(:, 1, :) - tmpsbs);
     
     tmpt2y = t(:, 1, :) - tmpu1y .* (t(:, 2, :) - t(:, 1, :));
     tmps2y = s(:, 1, :) - tmpu1y .* (s(:, 2, :) - s(:, 1, :));
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %       |1  1  0 ...   0  0|
     %       |0  0  1 ...   0  0|
     % R1 =  |0 -1  0 ...   1  0|
     %       |0  0 -1 ...   0  0|
     %       |0  0  0 ...  -1  1|
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     tmpwmy = AZF(w(:, 2, :)) .* dti ./ (Rk2(tmpzzy) .* tmpdty);
     tmpt3y = tmpt2y - tmpwmy .* (Rk2(t(:, 2, :)));
     tmps3y = tmps2y - tmpwmy .* (Rk2(s(:, 2, :)));
     
     tmpty = tmpt1y .* A1y + (tmpt2y .* B2y + B1y .* tmpt3y).*A2y;
     tmpsy = tmps1y .* A1y + (tmps2y .* B2y + B1y .* tmps3y).*A2y;
     
     uf(1:im, 1, 1:kbm1) = tmpty(:,1, 1:kbm1);
     vf(1:im, 1, 1:kbm1) = tmpsy(:,1, 1:kbm1); 
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     tmpfsm = repmat(fsm, 1, 1, kb);
     tmpfsm(:, :, kb) = 1;
     
     uf = uf .* tmpfsm;
     vf = vf .* tmpfsm;           
            
    return
    %
elseif(idx==5) 
    %     Vertical velocity boundary conditions:
%     tmpfsm = repmat(fsm, 1, 1, kb);
%     tmpfsm(:, :, kb) = 1;
%     w = w .* tmpfsm;   
    w=w.*fsm_3d;
    return
    %
elseif(idx==6)
    %
    %     q2 and q2l boundary conditions:
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     dx_3d = repmat(dx, 1, 1, kb); 
     tmpdx1 = dx_3d(im, :, :);
     tmpdx2 = dx_3d(imm1, :, :);

     
     tmpu1x = (2.e0 * u(im,:,:) * dti) ./ (tmpdx1 + tmpdx2);
     A1x = tmpu1x;
     A2x = tmpu1x;
     
%     tmp1=create_field((tmpu1x > 0.e0)
     A1x(double(tmpu1x > 0.e0)) = 0.e0;
     A1x(double(tmpu1x <= 0.e0)) = 1.e0;
     A2x(double(tmpu1x > 0.e0)) = 1.e0;
     A2x(double(tmpu1x <= 0.e0)) = 0.e0;
     
     tmpq21x = q2(im, :, :) - (tmpu1x .* (small - q2(im, :, :)));
     tmpq2l1x = q2l(im, :, :) - (tmpu1x .* (small - q2l(im, :, :)));
     
     tmpq22x = q2(im, :, :) - (tmpu1x .* (q2(im, :, :) - q2(imm1, :, :)));
     tmpq2l2x = q2l(im, :, :) - (tmpu1x .* (q2l(im, :, :) - q2l(imm1, :, :)));         
     
     uf(im, :, :) = A1x(:,:,:) .* tmpq21x(:,:,:)+ A2x(:,:,:) .* tmpq22x(:,:,:);
     vf(im, :, :) = A1x(:,:,:) .* tmpq2l1x(:,:,:)+ A2x(:,:,:) .* tmpq2l2x(:,:,:);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     
     tmpdx1 = dx_3d(1, :, :);
     tmpdx2 = dx_3d(2, :, :);
     
     tmpu1x = (2.e0 .* u(2, :, :) .* dti) ./ (tmpdx1 + tmpdx2);
     A1x = tmpu1x;
     A2x = tmpu1x;
     
     A1x(double(tmpu1x >= 0.e0)) = 1.e0;
     A1x(double(tmpu1x < 0.e0)) = 0.e0;
     A2x(double(tmpu1x >= 0.e0)) = 0.e0;
     A2x(double(tmpu1x < 0.e0)) = 1.e0;
     
     tmpq21x = q2(1, :, :) - tmpu1x .* (q2(1, :, :) - small);
     tmpq2l1x = q2l(1, :, :) - tmpu1x .* (q2l(1, :, :) - small);
     
     tmpq22x = q2(1, :, :) - tmpu1x .* (q2(2, :, :) - q2(1, :, :));
     tmpq2l2x = q2l(1, :, :) - tmpu1x .* (q2l(2, :, :) - q2l(1, :, :));
     
     uf(1, :, :) = A1x(1,:,:) .* tmpq21x(1,:,:) + A2x(1,:,:) .* tmpq22x(1,:,:);
     vf(1, :, :) = A1x(1,:,:) .* tmpq2l1x(1,:,:) + A2x(1,:,:) .* tmpq2l2x(1,:,:);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     dy_3d = repmat(dy, 1, 1, kb);
     tmpdy1 = dy_3d(:, jm, :);
     tmpdy2 = dy_3d(:, jmm1, :);
     
     tmpu1y = (2.e0 .* v(:, jm, :) .* dti) ./ (tmpdy1 + tmpdy2);
     A1y = tmpu1y;
     A2y = tmpu1y;
     
     A1y(double(tmpu1y > 0.e0)) = 0.e0;
     A1y(double(tmpu1y <= 0.e0)) = 1.e0;
     A2y(double(tmpu1y > 0.e0)) = 1.e0;
     A2y(double(tmpu1y <= 0.e0)) = 0.e0;
     
     tmpq21y = q2(:, jm, :) - tmpu1y .* (small - q2(:, jm, :));
     tmpq2l1y = q2l(:, jm, :) - tmpu1y .* (small - q2l(:, jm, :));
     
     tmpq22y = q2(:, jm, :) - tmpu1y .* (q2(:, jm, :) - q2(:, jmm1, :));
     tmpq2l2y = q2l(:, jm, :) - tmpu1y .* (q2l(:, jm, :) - q2l(:, jmm1, :));
     
     uf(:, jm, :) = tmpq21y(:,1,:) .* A1y(:,1,:) + tmpq22y(:,1,:) .* A2y(:,1,:);
     vf(:, jm, :) = tmpq2l1y(:,1,:) .* A1y(:,1,:) + tmpq2l2y(:,1,:) .* A2y(:,1,:);
       
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     tmpdy1 = dy_3d(:, 1, :);
     tmpdy2 = dy_3d(:, 2, :);
     
     tmpu1y = (2.e0 .* v(:, 2, :) .* dti) ./ (tmpdy1 + tmpdy2);
     A1y = tmpu1y;
     A2y = tmpu1y;
     
     A1y(double(tmpu1y >= 0.e0)) = 1.e0;
     A1y(double(tmpu1y < 0.e0)) = 0.e0;
     A2y(double(tmpu1y >= 0.e0)) = 0.e0;
     A2y(double(tmpu1y < 0.e0)) = 1.e0;
     
     tmpq21y = q2(:, 1, :) - tmpu1y .* (q2(:, 1, :) - small);
     tmpq2l1y = q2l(:, 1, :) - tmpu1y .* (q2l(:, 1, :) - small);
     
     tmpq22y = q2(:, 1, :) - tmpu1y .* (q2(:, 2, :) - q2(:, 1, :));
     tmpq2l2y = q2l(:, 1, :) - tmpu1y .* (q2l(:, 2, :) - q2l(:, 1, :));
     
     uf(:, 1, :) = tmpq21y(:,1,:) .* A1y(:,1,:) + tmpq22y(:,1,:) .* A2y(:,1,:);
     vf(:, 1, :) = tmpq2l1y(:,1,:) .* A1y(:,1,:) + tmpq2l2y(:,1,:) .* A2y(:,1,:);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     ssmall = 1.e-10;
     tmpfsm = repmat(fsm, 1, 1, kb);
     
     uf = uf .* tmpfsm + ssmall;
     vf = vf .* tmpfsm + ssmall;            

    return
end