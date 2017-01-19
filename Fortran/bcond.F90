subroutine new_bcond(idx, dti)
  use dm_op
  use grid
  use input
  use dm

  implicit none

  integer, intent(in)    :: idx
  real(kind=8), intent(in) :: dti
  integer                :: ierr
  type(Matrix)           :: tmpuvaf, local_array1, local_array2
  type(Matrix)           :: tmpuvf, tmpMASK
  type(Matrix)           :: tmpw
  type(Matrix)           :: tmp,tmpwm, tmp1, tmp2, tmp3, &
       tmp4, tmp5, tmp6, filter1, filter2, filter3, filter4
  type(Matrix)     :: Li_h, Li_el, Rj_h, Rj_el
  
  integer, allocatable :: seq_im(:), seq_jm(:)

  real(kind=8), allocatable ::  arr_im(:), arr_jm(:), arr_jmm1(:), arr_imm1(:), &
       arr_jmm2(:), arr_imm2(:)
  
  integer :: i, j
  
  allocate(arr_im(im), arr_jm(jm), seq_im(im), seq_jm(jm))
  allocate(arr_imm1(imm1), arr_jmm1(jmm1), arr_imm2(imm2), arr_jmm2(jmm2))
  
  seq_im = (/(0, i=0,im-1)/)
  seq_jm = (/(0, j=0,jm-1)/)

  if (idx == 1) then

     ! elf(1, :) = elf(2, :);
     call dm_getvalues(elf, (/1/), seq_jm, (/0/), arr_jm, ierr)
     call dm_setvalues(elf, (/0/), seq_jm, (/0/), arr_jm, ierr)

     ! elf(im, :) = elf(imm1, :);
     call dm_getvalues(elf, (/imm2/), seq_jm, (/0/), arr_jm, ierr)
     call dm_setvalues(elf, (/imm1/), seq_jm, (/0/), arr_jm, ierr)

     ! elf(:, 1) = elf(:, 2);
     call dm_getvalues(elf, seq_im, (/1/), (/0/), arr_im, ierr)
     call dm_setvalues(elf, seq_im, (/0/), (/0/), arr_im, ierr)

     ! elf(:, jm) = elf(:, jmm1);
     call dm_getvalues(elf, seq_im, (/jmm2/), (/0/), arr_im, ierr)
     call dm_setvalues(elf, seq_im, (/jmm1/), (/0/), arr_im, ierr)

     ! elf = elf .* fsm; 
     elf = elf .em. fsm

     ! elf = elf .em. REV_MASK_X1 .em. REV_MASK_X2 +
     !       shift(elf, 1, 1) .em. MASK_X1 + shift(elf, 1, -1) .em. MASK_X2
     ! elf = elf .em. REV_MASK_Y1 .em. REV_MASK_Y2 +
     !       shift(elf, 2, 1) .em. MASK_Y1 + shift(elf, 2, -1) .em. MASK_Y2
     ! elf = elf .em. fsm

  else if (idx == 2) then

     ! tmph = h;
     ! tmpel = el;
     ! tmph(im, :) = tmph(imm1, :);
     ! tmpel(im, :) = tmpel(imm1, :);
     ! uab1 = repmat(uabe, im, 1);
     ! el1 = repmat(ele, im, 1);
     ! tmpuaf = uab1 + rfe .* sqrt(grav ./(tmph)) .* (tmpel - el1);
     ! tmpuaf = ramp .* tmpuaf;
     ! uaf(im, 2:jmm1) = tmpuaf(im, 2:jmm1);
     !
     ! 
     ! tmph(1, :) = tmph(2, :);
     ! tmpel(1, :) = tmpel(2, :);
     ! uab2 = repmat(uabw, im, 1);
     ! el2 = repmat(elw, im, 1);
     ! tmpuaf = uab2 - rfw .* sqrt(grav ./(tmph)) .* (tmpel - el2);
     ! tmpuaf = ramp .* tmpuaf;
     ! uaf(1, 2:jmm1) = tmpuaf(1, 2:jmm1);
     ! uaf(2, 2:jmm1) = tmpuaf(1, 2:jmm1);
     ! 
     ! vaf(im, 2:jmm1) = 0.e0;
     ! vaf(1, 2:jmm1) = 0.e0;
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! tmph = h;
     ! tmpel = el;
     ! tmph(:, jm) = tmph(:, jmm1);
     ! tmpel(:, jm) = tmpel(:, jmm1);
     ! vab1 = repmat(vabn, jm, 1);
     ! vab1 = vab1';
     ! el3 = repmat(eln, jm, 1);
     ! el3 = el3';
     ! tmpvaf = vab1 + rfn .* sqrt(grav./(tmph)).*(tmpel - el3);
     ! tmpvaf = ramp .* tmpvaf;
     ! vaf(2:imm1, jm) = tmpvaf(2:imm1, jm);
     ! 
     ! tmph(:, 1) = tmph(:, 2);
     ! tmpel(:, 1) = tmpel(:, 2);
     ! vab2 = repmat(vabs, jm, 1);
     ! vab2 = vab2';
     ! el4 = repmat(els, jm, 1);
     ! el4 = el4';
     ! tmpvaf = vab2 - rfs .* sqrt(grav./(tmph)).*(tmpel - el4);
     ! tmpvaf = ramp .* tmpvaf;
     ! vaf(2:imm1, 1) = tmpvaf(2:imm1, 1);
     ! vaf(2:imm1, 2) = tmpvaf(2:imm1, 1);
     ! 
     ! uaf(2:imm1, jm) = 0.e0;
     ! uaf(2:imm1, 1) = 0.e0;
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! uaf(:,:) = uaf(:,:) .* dum(:,:);
     ! vaf(:,:) = vaf(:,:) .* dvm(:,:);

     
     Li_el = COPY_EDGE_X(el)
     Li_h = COPY_EDGE_X(h)
     
     tmpuvaf = dm_rep(uabe, im, 1, 1) + &
          rfe * dm_sqrt(grav * dm_pow(Li_h, -1)) .em. (Li_el - dm_rep(ele, im, 1, 1))

     tmpuvaf = ramp * tmpuvaf

     call dm_getvalues(tmpuvaf, (/imm1/), (/(j, j=1,jmm2)/), (/0/), arr_jmm2, ierr)
     call dm_setvalues(uaf, (/imm1/), (/(j, j=1,jmm2)/), (/0/), arr_jmm2, ierr)

     tmpuvaf = dm_rep(uabw, im, 1, 1) - &
          rfw * dm_sqrt(grav * dm_pow(Li_h, -1)) .em. (Li_el - dm_rep(elw, im, 1, 1))

     tmpuvaf = ramp * tmpuvaf

     call dm_getvalues(tmpuvaf, (/0/), (/(j, j=1,jmm2)/), (/0/), arr_jmm2, ierr)
     call dm_setvalues(uaf, (/0/), (/(j, j=1,jmm2)/),(/0/),  arr_jmm2, ierr)
     call dm_setvalues(uaf, (/1/), (/(j, j=1,jmm2)/), (/0/), arr_jmm2, ierr)

     call dm_setvalues(vaf, (/imm1/), (/(j, j=1,jmm2)/), (/0/), arr_jmm2, ierr)
     call dm_setvalues(vaf, (/0/), (/(j, j=1,jmm2)/),(/0/), arr_jmm2, ierr)

     Rj_h = COPY_EDGE_Y(h)
     Rj_el = COPY_EDGE_Y(el)
     tmpuvaf = dm_trans(dm_rep(vabn, jm, 1, 1)) + &
          rfn * dm_sqrt(grav * dm_pow(Rj_h, -1)) .em. (Rj_el - &
          dm_trans(dm_rep(eln, jm, 1, 1)))

     tmpuvaf = ramp * tmpuvaf

     call dm_getvalues(tmpuvaf, (/(i, i=1,imm2)/), (/jmm1/), (/0/), arr_imm2, ierr)
     call dm_setvalues(vaf, (/(i, i=1,imm2)/), (/jmm1/),(/0/), arr_imm2, ierr)

     tmpuvaf = dm_trans(dm_rep(vabs, jm, 1, 1)) - &
          rfs * dm_sqrt(grav * dm_pow(Rj_h, -1)) .em. (Rj_el - &
          dm_trans(dm_rep(els, jm, 1, 1)))

     tmpuvaf = ramp * tmpuvaf

     call dm_getvalues(tmpuvaf, (/(i, i=1,imm2)/), (/0/), (/0/), arr_imm2, ierr)
     call dm_setvalues(vaf, (/(i, i=1,imm2)/), (/0/),(/0/), arr_imm2, ierr)
     call dm_setvalues(vaf, (/(i, i=1,imm2)/), (/1/),(/0/), arr_imm2, ierr)

     call dm_setvalues(uaf, (/(i, i=1,imm2)/), (/jmm1/),(/0/), arr_imm2, ierr)
     call dm_setvalues(uaf, (/(i, i=1,imm2)/), (/0/),(/0/), arr_imm2, ierr)

     uaf = uaf .em. dum
     vaf = vaf .em. dvm

     call dm_destroy(tmpuvaf, ierr)
     call dm_destroy(Li_h, ierr)
     call dm_destroy(Li_el, ierr)
     call dm_destroy(Rj_h, ierr)
     call dm_destroy(Rj_el, ierr)
  else if (idx == 3) then

     ! h_3d = repmat(h, 1, 1, kb);
     ! tmpu1 = u(imm1, :, :);
     ! tmpu2 = u(im, :, :);
     ! tmpga = sqrt(h_3d(im, :, :) ./ hmax);
     ! 
     ! tmpuf = tmpga .* (AXF1(AXB1(tmpu1))) + (1.e0 - tmpga) .* (AXF1(AXB1(tmpu2)));
     ! uf(im, 2:jmm1, 1:kbm1) = tmpuf(1, 2:jmm1, 1:kbm1);
     ! vf(im, 2:jmm1, 1:kbm1) = 0.e0;
     ! 
     ! tmpu1 = u(3, :, :);
     ! tmpu2 = u(2, :, :);
     ! tmpga = sqrt(h_3d(1, :, :) ./ hmax);
     ! 
     ! tmpuf = tmpga .* (AXF1(AXB1(tmpu1))) + (1.e0 - tmpga) .* (AXF1(AXB1(tmpu2)));
     ! uf(2, 2:jmm1, 1:kbm1) = tmpuf(1, 2:jmm1, 1:kbm1);
     ! uf(1, 2:jmm1, 1:kbm1) = tmpuf(1, 2:jmm1, 1:kbm1);
     ! vf(1, 2:jmm1, 1:kbm1) = 0.e0;
     ! 
     !%%          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

     ! tmpv1 = v(:, jmm1, :);
     ! tmpv2 = v(:, jm, :);
     ! tmpga = sqrt(h_3d(:, jm, :) ./ hmax);
     ! 
     ! tmpvf = tmpga .* (AXF2(AXB2(tmpv1))) + (1.e0 - tmpga) .* (AXF2(AXB2(tmpv2)));
     ! vf(2:imm1, jm, 1:kbm1) = tmpvf(2:imm1,1,1:kbm1);
     ! uf(2:imm1, jm, 1:kbm1) = 0.e0;
     ! 
     ! tmpv1 = v(:, 3, :);
     ! tmpv2 = v(:, 2, :);
     ! tmpga = sqrt(h_3d(:, 1, :) ./ hmax);
     ! 
     ! tmpvf = tmpga .* (AXF2(AXB2(tmpv1))) + (1.e0 - tmpga) .* (AXF2(AXB2(tmpv2)));
     ! vf(2:imm1, 2, 1:kbm1) = tmpvf(2:imm1,1,1:kbm1);
     ! vf(2:imm1, 1, 1:kbm1) = tmpvf(2:imm1,1,1:kbm1);
     ! uf(2:imm1, 1, 1:kbm1) = 0.e0;

     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                tmpdum=repmat(dum,1,1,kb);
     ! tmpdvm=repmat(dvm,1,1,kb);
     !
     ! tmpdum(:, :, kb) = 1.e0;
     ! tmpdvm(:, :, kb) = 1.e0;
     !
     ! uf = uf .* tmpdum;
     ! vf = vf .* tmpdvm;

     h_3d = dm_rep(h, 1, 1, kb)

     tmpuvf = dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_X2)) .em. &
          (AXF1(AXB1(shift(u, 1, -1) .em. MASK_X2))) + &
          (1.e0 - dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_X2))) .em. &
          (AXF1(AXB1(u .em. MASK_X2)))

     tmpMASK = MASK_Z2
     call dm_setvalues(tmpMASK, (/imm1/), (/0/), (/kbm1/), (/0.e0/), ierr)
     call dm_setvalues(tmpMASK, (/imm1/), (/jmm1/), (/kbm1/), (/0.e0/), ierr)
     call dm_setvalues(tmpMASK, (/0/), (/0/), (/kbm1/), (/0.e0/), ierr)
     call dm_setvalues(tmpMASK, (/0/), (/jmm1/), (/kbm1/), (/0.e0/), ierr)

     uf   = uf .em. REV_MASK_X2 + uf .em. MASK_X2 .em. MASK_Y1 +&
          uf .em. MASK_X2 .em. MASK_Y2 + &
          uf .em. MASK_X2 .em. tmpMASK + &
          tmpuvf .em. REV_MASK_Y1 .em. REV_MASK_Y2 .em. REV_MASK_Z2

     tmpuvf = dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_X1)) .em. &
          (AXF1(AXB1(shift(shift(u, 1, 1), 1, 1) .em. MASK_X1))) +&
          (1.e0 - dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_X1))) .em. &
          (AXF1(AXB1(shift(u, 1, 1) .em. MASK_X1)))

     uf   = uf .em. REV_MASK_X1 + uf .em. MASK_X1 .em. MASK_Y1 + &
          uf .em. MASK_X1 .em. MASK_Y2 + &
          uf .em. MASK_X1 .em. tmpMASK + tmpuvf .em. REV_MASK_Y1 .em. &
          REV_MASK_Y2 .em. REV_MASK_Z2

     uf   = uf .em. MASK_X1 + &
          shift(shift(shift(shift(uf, 1, 1), 1, 1), 1, -1), 1, -1) + &
          shift((shift(uf, 1, 1) .em. MASK_X1 .em. MASK_Y1 + &
          shift(uf, 1, 1) .em. MASK_X1 .em. MASK_Y2 + &
          shift(uf, 1, 1) .em. MASK_X1 .em. tmpMASK), 1, -1) + &
          shift((tmpuvf .em. REV_MASK_Y1 .em. REV_MASK_Y2 .em. REV_MASK_Z2), 1, -1)

     vf = vf .em. REV_MASK_X1 .em. REV_MASK_X2 + vf .em. &
          (MASK_X1 + MASK_X2) .em. (MASK_Y1 + MASK_Y2) + &
          vf .em. (MASK_X1 + MASK_X2) .em. tmpMASK

     tmpuvf= dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_Y2)) .em. &
          (AXF2(AXB2(shift(v, 1, -1) .em. MASK_Y2))) + &
          (1.e0 - dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_Y2))) .em. &
          (AXF2(AXB2(v .em. MASK_Y2)))

     ! tmpMASK = MASK_Z2
     ! call dm_setvalue(tmpMASK, 1, jm, kb, 0.e0, ierr)
     ! call dm_setvalue(tmpMASK, im, jm, kb, 0.e0, ierr)

     vf   = vf .em. REV_MASK_Y2 + vf .em. MASK_Y2 .em. MASK_X1 +&
          vf .em. MASK_Y2 .em. MASK_X2 + &
          vf .em. MASK_Y2 .em. tmpMASK + &
          tmpuvf .em. REV_MASK_X1 .em. REV_MASK_X2 .em. REV_MASK_Z2

     tmpuvf= dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_Y1)) .em. &
          (AXF2(AXB2(shift(shift(v, 2, 1), 2, 1) .em. MASK_Y1))) +&
          (1.e0 - dm_sqrt((1.e0/hmax) * (h_3d .em. MASK_Y1))) .em. &
          (AXF2(AXB2(shift(v, 2, 1) .em. MASK_Y1)))

     ! call dm_setvalue(tmpMASK, 1, 1, kb, 0.e0, ierr)
     ! call dm_setvalue(tmpMASK, im, 1, kb, 0.e0, ierr)

     vf   = vf .em. REV_MASK_Y1 + vf .em. MASK_Y1 .em. MASK_X1 + &
          vf .em. MASK_Y1 .em. MASK_X2 + &
          vf .em. MASK_Y1 .em. tmpMASK + &
          tmpuvf .em. REV_MASK_X1 .em. REV_MASK_X2 .em. REV_MASK_Z2

     vf   = vf .em. MASK_Y1 + &
          shift(shift(shift(shift(vf, 2, 1), 2, 1), 2, -1), 2, -1) + &
          shift((shift(vf, 2, 1) .em. MASK_Y1 .em. MASK_X1 + &
          shift(vf, 2, 1) .em. MASK_Y1 .em. MASK_X2 + &
          shift(vf, 2, 1) .em. MASK_Y1 .em. tmpMASK), 2, -1) + &
          shift((tmpuvf .em. REV_MASK_X1 .em. REV_MASK_X2 .em. REV_MASK_Z2), 2, -1)

     uf = uf .em. REV_MASK_Y1 .em. REV_MASK_Y2 + &
          uf .em. (MASK_Y1 + MASK_Y2) .em. (MASK_X1 + MASK_X2) + &
          uf .em. (MASK_Y1 + MASK_Y2) .em. tmpMASK

     uf = uf .em. MASK_Z2 + uf .em. dm_rep(dum, 1, 1, kb) .em. REV_MASK_Z2
     vf = vf .em. MASK_Z2 + vf .em. dm_rep(dvm, 1, 1, kb) .em. REV_MASK_Z2

     call dm_destroy(tmpuvf, ierr)
     call dm_destroy(tmpMASK, ierr)

  else if (idx == 4) then

     ! dx_3d = repmat(dx, 1, 1, kb);
     ! dt_3d = repmat(dt, 1, 1, kb);
     ! tmpdtx = dt_3d(imm1, :, :);
     ! tmpdx1 = dx_3d(im, :, :);
     ! tmpdx2 = dx_3d(imm1, :, :);  
     ! tmpzzx = zeros(1, jm, kb);
     ! tmpzzx(1,:,:) = repmat(zz, jm, 1);
     ! 
     ! B1x = zeros(1, jm, kb);
     ! B2x = zeros(1, jm ,kb);
     ! 
     ! tmpu1x = (2.e0 .* u(im, :, :) .* dti) ./ (tmpdx1 + tmpdx2);
     ! A1x = tmpu1x;
     ! A2x = tmpu1x;
     ! 
     ! A1x(tmpu1x > 0.e0) = 0.e0;
     ! A1x(tmpu1x <= 0.e0) = 1.e0;
     ! A2x(tmpu1x > 0.e0) = 1.e0;
     ! A2x(tmpu1x <= 0.e0) = 0.e0;
     ! B1x(:,:,2:kb-2) = 1.e0;
     ! B2x(B1x < .9e0) = 1.e0;

     ! tmptbe(1,:,:) = tbe;
     ! tmpsbe(1,:,:) = sbe;
     ! 
     ! tmpt1x = t(im, :, :) - tmpu1x .* (tmptbe - t(im, :, :));
     ! tmps1x = s(im, :, :) - tmpu1x .* (tmpsbe - s(im, :, :));
     ! 
     ! tmpt2x = t(im, :, :) - tmpu1x .* (t(im, :, :) - t(imm1, :, :));
     ! tmps2x = s(im, :, :) - tmpu1x .* (s(im, :, :) - s(imm1, :, :));
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! %       |1  1  0 ...   0  0|
     ! %       |0  0  1 ...   0  0|
     ! % R1 =  |0 -1  0 ...   1  0|
     ! %       |0  0 -1 ...   0  0|
     ! %       |0  0  0 ...  -1  1|
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! tmpwmx = AZF1(w(imm1, :, :)) .* dti ./ (Rk1(tmpzzx) .* tmpdtx);
     ! tmpt3x = tmpt2x - tmpwmx .*  (Rk1(t(imm1, :, :)));
     ! tmps3x = tmps2x - tmpwmx .* (Rk1(s(imm1, :, :)));
     ! 
     ! tmptx = tmpt1x .* A1x + (tmpt2x .* B2x + B1x .* tmpt3x).*A2x;
     ! tmpsx = tmps1x .* A1x + (tmps2x .* B2x + B1x .* tmps3x).*A2x;
     ! 
     ! uf(im, 1:jm, 1:kbm1) = tmptx(1, :, 1:kbm1);
     ! vf(im, 1:jm, 1:kbm1) = tmpsx(1, :, 1:kbm1);

     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! tmpdtx = dt_3d(2, :, :);   
     ! tmpdx1 = dx_3d(1, :, :);
     ! tmpdx2 = dx_3d(2, :, :);
     ! 
     ! tmpu1x = (2.e0 .* u(2, :, :) .* dti) ./ (tmpdx1 + tmpdx2);
     ! A1x = tmpu1x;
     ! A2x = tmpu1x;
     ! 
     ! A1x(tmpu1x >= 0.e0) = 1.e0;
     ! A1x(tmpu1x < 0.e0) = 0.e0;
     ! A2x(tmpu1x >= 0.e0) = 0.e0;
     ! A2x(tmpu1x < 0.e0) = 1.e0;
     ! 
     ! tmptbw(1,:,:) = tbw;
     ! tmpsbw(1,:,:) = sbw;
     ! 
     ! tmpt1x = t(1, :, :) - tmpu1x .* (t(1, :, :) - tmptbw);
     ! tmps1x = s(1, :, :) - tmpu1x .* (s(1, :, :) - tmpsbw);
     ! 
     ! tmpt2x = t(1, :, :) - tmpu1x .* (t(2, :, :) - t(1, :, :));
     ! tmps2x = s(1, :, :) - tmpu1x .* (s(2, :, :) - s(1, :, :));
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! %       |1  1  0 ...   0  0|
     ! %       |0  0  1 ...   0  0|
     ! % R1 =  |0 -1  0 ...   1  0|
     ! %       |0  0 -1 ...   0  0|
     ! %       |0  0  0 ...  -1  1|
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! tmpwmx = AZF1(w(2, :, :)) .* dti ./ (Rk1(tmpzzx) .* tmpdtx);
     ! tmpt3x = tmpt2x - tmpwmx .* (Rk1(t(2, :, :)));
     ! tmps3x = tmps2x - tmpwmx .* (Rk1(s(2, :, :)));
     ! 
     ! tmptx = tmpt1x .* A1x + (tmpt2x .* B2x + B1x .* tmpt3x).*A2x;
     ! tmpsx = tmps1x .* A1x + (tmps2x .* B2x + B1x .* tmps3x).*A2x;
     ! 
     ! uf(1, 1:jm, 1:kbm1) = tmptx(1, :, 1:kbm1);
     ! vf(1, 1:jm, 1:kbm1) = tmpsx(1, :, 1:kbm1);
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! dy_3d = repmat(dy, 1, 1, kb);
     ! tmpdty = dt_3d(:, jmm1, :); 
     ! tmpdy1 = dy_3d(:, jm, :);
     ! tmpdy2 = dy_3d(:, jmm1, :);
     ! tmpzzy(:,1,:) = repmat(zz, im, 1);

     ! B1y = zeros(im, 1, kb);
     ! B2y = zeros(im, 1, kb);

     ! tmpu1y = (2.e0 .* v(:, jm, :) .* dti) ./ (tmpdy1 + tmpdy2);
     ! A1y = tmpu1y;
     ! A2y = tmpu1y;
     ! 
     ! A1y(tmpu1y > 0.e0) = 0.e0;
     ! A1y(tmpu1y <= 0.e0) = 1.e0;
     ! A2y(tmpu1y > 0.e0) = 1.e0;
     ! A2y(tmpu1y <= 0.e0) = 0.e0;
     ! B1y(:,:, 2:kb-2) = 1;
     ! B2y(B1y < .9e0) = 1;

     ! tmptbn(:,1,:) = tbn;
     ! tmpsbn(:,1,:) = sbn;
     ! 
     ! tmpt1y = t(:, jm, :) - tmpu1y .* (tmptbn - t(:, jm, :));
     ! tmps1y = s(:, jm, :) - tmpu1y .* (tmpsbn - s(:, jm, :));
     ! 
     ! tmpt2y = t(:, jm, :) - tmpu1y .* (t(:, jm, :) - t(:, jmm1, :));
     ! tmps2y = s(:, jm, :) - tmpu1y .* (s(:, jm, :) - s(:, jmm1, :));
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! %       |1  1  0 ...   0  0|
     ! %       |0  0  1 ...   0  0|
     ! % R1 =  |0 -1  0 ...   1  0|
     ! %       |0  0 -1 ...   0  0|
     ! %       |0  0  0 ...  -1  1|
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! tmpwmy = AZF(w(:, jmm1, :)) .* dti ./ (Rk2(tmpzzy) .* tmpdty);
     ! tmpt3y = tmpt2y - tmpwmy .* (Rk2(t(:, jmm1, :)));
     ! tmps3y = tmps2y - tmpwmy .* (Rk2(s(:, jmm1, :)));
     ! 
     ! tmpty = tmpt1y .* A1y + (tmpt2y .* B2y + B1y .* tmpt3y).*A2y;
     ! tmpsy = tmps1y .* A1y + (tmps2y .* B2y + B1y .* tmps3y).*A2y;

     ! uf(1:im, jm, 1:kbm1) = tmpty(:,1, 1:kbm1);
     ! vf(1:im, jm, 1:kbm1) = tmpsy(:,1, 1:kbm1);       
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ! tmpdty = dt_3d(:, 2, :);     
     ! tmpdy1 = dy_3d(:, 1, :);
     ! tmpdy2 = dy_3d(:, 2, :);
     ! 
     ! tmpu1y = (2.e0 .* v(:, 2, :) .* dti) ./ (tmpdy1 + tmpdy2);
     ! A1y = tmpu1y;
     ! A2y = tmpu1y;
     ! 
     ! A1y(tmpu1y >= 0.e0) = 1.e0;
     ! A1y(tmpu1y < 0.e0) = 0.e0;
     ! A2y(tmpu1y >= 0.e0) = 0.e0;
     ! A2y(tmpu1y < 0.e0) = 1.e0;
     ! 
     ! tmptbs(:,1,:) = tbs;
     ! tmpsbs(:,1,:) = sbs;
     ! 
     ! tmpt1y = t(:, 1, :) - tmpu1y .* (t(:, 1, :) - tmptbs);
     ! tmps1y = s(:, 1, :) - tmpu1y .* (s(:, 1, :) - tmpsbs);
     ! 
     ! tmpt2y = t(:, 1, :) - tmpu1y .* (t(:, 2, :) - t(:, 1, :));
     ! tmps2y = s(:, 1, :) - tmpu1y .* (s(:, 2, :) - s(:, 1, :));
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! %       |1  1  0 ...   0  0|
     ! %       |0  0  1 ...   0  0|
     ! % R1 =  |0 -1  0 ...   1  0|
     ! %       |0  0 -1 ...   0  0|
     ! %       |0  0  0 ...  -1  1|
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! tmpwmy = AZF(w(:, 2, :)) .* dti ./ (Rk2(tmpzzy) .* tmpdty);
     ! tmpt3y = tmpt2y - tmpwmy .* (Rk2(t(:, 2, :)));
     ! tmps3y = tmps2y - tmpwmy .* (Rk2(s(:, 2, :)));
     ! 
     ! tmpty = tmpt1y .* A1y + (tmpt2y .* B2y + B1y .* tmpt3y).*A2y;
     ! tmpsy = tmps1y .* A1y + (tmps2y .* B2y + B1y .* tmps3y).*A2y;
     ! 
     ! uf(1:im, 1, 1:kbm1) = tmpty(:,1, 1:kbm1);
     ! vf(1:im, 1, 1:kbm1) = tmpsy(:,1, 1:kbm1); 
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! tmpfsm = repmat(fsm, 1, 1, kb);
     ! tmpfsm(:, :, kb) = 1;
     ! 
     ! uf = uf .* tmpfsm;
     ! vf = vf .* tmpfsm;           



     dx_3d = dm_rep(dx, 1, 1, kb)
     dt_3d = dm_rep(dt, 1, 1, kb)
     dy_3d = dm_rep(dy, 1, 1, kb)

     tmp = (2.e0 * dti * (u .em. MASK_X2)) .ed. &
          (dx_3d .em. MASK_X2 + shift(dx_3d, 1, -1) .em. MASK_X2)

     filter1 = (tmp <= 0.e0)
     filter2 = (tmp > 0.e0)
     filter3 = shift(shift(shift(shift(dm_ones(im, jm, kb), 3, -1), 3, -1), 3, 1), 3, 1) .em. REV_MASK_Z1
     filter4 = (filter3 < 0.9e0)

     filter1 = filter1 .em. MASK_X2 
     filter2 = filter2 .em. MASK_X2 
     filter3 = filter3 .em. MASK_X2 
     filter4 = filter4 .em. MASK_X2

     tmp1 = t .em. MASK_X2 - tmp .em. (dm_rep(tbe, im, 1, 1) .em. &
          MASK_X2 - t .em. MASK_X2)
     tmp2 = s .em. MASK_X2 - tmp .em. (dm_rep(sbe, im, 1, 1) .em. &
          MASK_X2 - s .em. MASK_X2)

     tmp3 = t .em. MASK_X2 - tmp .em. &
          (t .em. MASK_X2 - shift(t, 1, -1) .em. MASK_X2)
     tmp4 = s .em. MASK_X2 - tmp .em. &
          (s. em. MASK_X2 - shift(s, 1, -1) .em. MASK_X2)

     tmpwm= dti * AZF1(shift(w, 1, -1) .em. MASK_X2) .ed. &
          (Rk1(zz_3d .em. MASK_X2) .em. (shift(dt_3d, 1, -1) .em. MASK_X2))

     tmp5 = tmp3 - tmpwm .em. (Rk1(shift(t, 1, -1) .em. MASK_X2))
     tmp6 = tmp4 - tmpwm .em. (Rk1(shift(s, 1, -1) .em. MASK_X2))

     uf = uf .em. REV_MASK_X2 + uf .em. MASK_X2 .em. MASK_Z2 + &
          (tmp1 .em. filter1 + (tmp3 .em. filter4 + tmp5 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2
     vf = vf .em. REV_MASK_X2 + vf .em. MASK_X2 .em. MASK_Z2 + &
          (tmp2 .em. filter1 + (tmp4 .em. filter4 + tmp6 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2

     tmp = (2.e0 * dti * (shift(u, 1, 1) .em. MASK_X1)) .ed. &
          (dx_3d .em. MASK_X1 + shift(dx_3d, 1, 1) .em. MASK_X1)

     filter1 = (tmp >= 0.e0)
     filter2 = (tmp < 0.e0)
     filter3 = shift(shift(shift(shift(dm_ones(im, jm, kb), 3, -1), 3, -1), 3, 1), 3, 1) .em. REV_MASK_Z1
     filter4 = (filter3 < 0.9e0)

     filter1 = filter1 .em. MASK_X1 
     filter2 = filter2 .em. MASK_X1 
     filter3 = filter3 .em. MASK_X1 
     filter4 = filter4 .em. MASK_X1


     tmp1 = t .em. MASK_X1 - tmp .em. &
          (t .em. MASK_X1 - dm_rep(tbw, im, 1, 1) .em. MASK_X1)
     tmp2 = s .em. MASK_X1 - tmp .em. &
          (s .em. MASK_X1 - dm_rep(sbw, im, 1, 1) .em. MASK_X1)

     tmp3 = t .em. MASK_X1 - tmp .em. (shift(t, 1, 1) .em. MASK_X1 - t .em. MASK_X1)
     tmp4 = s .em. MASK_X1 - tmp .em. (shift(s, 1, 1). em. MASK_X1 - s .em. MASK_X1)

     tmpwm= dti * AZF1(shift(w, 1, 1) .em. MASK_X1) .ed. &
          (Rk1(zz_3d .em. MASK_X1) .em. (shift(dt_3d, 1, 1) .em. MASK_X1))

     tmp5 = tmp3 - tmpwm .em. (Rk1(shift(t, 1, 1) .em. MASK_X1))
     tmp6 = tmp4 - tmpwm .em. (Rk1(shift(s, 1, 1) .em. MASK_X1))

     uf = uf .em. REV_MASK_X1 + uf .em. MASK_X1 .em. MASK_Z2 + &
          (tmp1 .em. filter1 + (tmp3 .em. filter4 + tmp5 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2
     vf = vf .em. REV_MASK_X1 + vf .em. MASK_X1 .em. MASK_Z2 + &
          (tmp2 .em. filter1 + (tmp4 .em. filter4 + tmp6 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2

     tmp = (2.e0 * dti * (v .em. MASK_Y2)) .ed. &
          (dy_3d .em. MASK_Y2 + shift(dy_3d, 2, -1) .em. MASK_Y2)

     filter1 = (tmp <= 0.e0)
     filter2 = (tmp > 0.e0)
     filter3 = shift(shift(shift(shift(dm_ones(im, jm, kb), 3, -1), 3, -1), 3, 1), 3, 1) .em. REV_MASK_Z1
     filter4 = (filter3 < 0.9e0)

     filter1 = filter1 .em. MASK_Y2 
     filter2 = filter2 .em. MASK_Y2 
     filter3 = filter3 .em. MASK_Y2 
     filter4 = filter4 .em. MASK_Y2

     tmp1 = t .em. MASK_Y2 - tmp .em. &
          (dm_rep(tbn, 1, jm, 1) .em. MASK_Y2 - t .em. MASK_Y2)
     tmp2 = s .em. MASK_Y2 - tmp .em. &
          (dm_rep(sbn, 1, jm, 1) .em. MASK_Y2 - s .em. MASK_Y2)

     tmp3 = t .em. MASK_Y2 - tmp .em. (t .em. MASK_Y2 - shift(t, 2, -1) .em. MASK_Y2)
     tmp4 = s .em. MASK_Y2 - tmp .em. (s. em. MASK_Y2 - shift(s, 2, -1) .em. MASK_Y2)

     tmpwm= dti * AZF(shift(w, 2, -1) .em. MASK_Y2) .ed. &
          (Rk2(zz_3d .em. MASK_Y2) .em. (shift(dt_3d, 2, -1) .em. MASK_Y2))

     tmp5 = tmp3 - tmpwm .em. (Rk2(shift(t, 2, -1) .em. MASK_Y2))
     tmp6 = tmp4 - tmpwm .em. (Rk2(shift(s, 2, -1) .em. MASK_Y2))

     uf = uf .em. REV_MASK_Y2 + uf .em. MASK_Y2 .em. MASK_Z2 + &
          (tmp1 .em. filter1 + (tmp3 .em. filter4 + tmp5 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2
     vf = vf .em. REV_MASK_Y2 + vf .em. MASK_Y2 .em. MASK_Z2 + &
          (tmp2 .em. filter1 + (tmp4 .em. filter4 + tmp6 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2

     tmp = (2.e0 * dti * (shift(v, 2, 1) .em. MASK_Y1)) .ed. &
          (dy_3d .em. MASK_Y1 + shift(dy_3d, 2, 1) .em. MASK_Y1)

     filter1 = (tmp >= 0.e0)
     filter2 = (tmp < 0.e0)
     filter3 = shift(shift(shift(shift(dm_ones(im, jm, kb), 3, -1), 3, -1), 3, 1), 3, 1) .em. REV_MASK_Z1
     filter4 = (filter3 < 0.9e0)

     filter1 = filter1 .em. MASK_Y1 
     filter2 = filter2 .em. MASK_Y1 
     filter3 = filter3 .em. MASK_Y1 
     filter4 = filter4 .em. MASK_Y1

     tmp1 = t .em. MASK_Y1 - tmp .em. &
          (t .em. MASK_Y1 - dm_rep(tbs, 1, jm, 1) .em. MASK_Y1)
     tmp2 = s .em. MASK_Y1 - tmp .em. &
          (s .em. MASK_Y1 - dm_rep(sbs, 1, jm, 1) .em. MASK_Y1)

     tmp3 = t .em. MASK_Y1 - tmp .em. (shift(t, 2, 1) .em. MASK_Y1 - t .em. MASK_Y1)
     tmp4 = s .em. MASK_Y1 - tmp .em. (shift(s, 2, 1). em. MASK_Y1 - s .em. MASK_Y1)

     tmpwm= dti * AZF(shift(w, 2, 1) .em. MASK_Y1) .ed. &
          (Rk2(zz_3d .em. MASK_Y1) .em. (shift(dt_3d, 2, 1) .em. MASK_Y1))

     tmp5 = tmp3 - tmpwm .em. (Rk2(shift(t, 2, 1) .em. MASK_Y1))
     tmp6 = tmp4 - tmpwm .em. (Rk2(shift(s, 2, 1) .em. MASK_Y1))

     uf = uf .em. REV_MASK_Y1 + uf .em. MASK_Y1 .em. MASK_Z2 + &
          (tmp1 .em. filter1 + (tmp3 .em. filter4 + tmp5 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2
     vf = vf .em. REV_MASK_Y1 + vf .em. MASK_Y1 .em. MASK_Z2 + &
          (tmp2 .em. filter1 + (tmp4 .em. filter4 + tmp6 .em.filter3) .em. filter2) &
          .em. REV_MASK_Z2

     uf = uf .em. MASK_Z2 + uf .em. dm_rep(fsm, 1, 1, kb) .em. REV_MASK_Z2
     vf = vf .em. MASK_Z2 + vf .em. dm_rep(fsm, 1, 1, kb) .em. REV_MASK_Z2

     call dm_destroy(tmp, ierr)
     call dm_destroy(tmpwm, ierr)
     call dm_destroy(tmp1, ierr)
     call dm_destroy(tmp2, ierr)
     call dm_destroy(tmp3, ierr)
     call dm_destroy(tmp4, ierr)
     call dm_destroy(tmp5, ierr)
     call dm_destroy(tmp6, ierr)
     call dm_destroy(dx_3d, ierr)
     call dm_destroy(dy_3d, ierr)
     call dm_destroy(dt_3d, ierr)
     call dm_destroy(filter1, ierr)
     call dm_destroy(filter2, ierr)
     call dm_destroy(filter3, ierr)
     call dm_destroy(filter4, ierr)

  else if (idx == 5) then

     ! tmpfsm = repmat(fsm, 1, 1, kb);
     ! tmpfsm(:, :, kb) = 1;
     ! w = w .* tmpfsm;   

     tmpw = w .em. dm_rep(fsm, 1, 1, kb)

     w = tmpw .em. REV_MASK_Z2 + w .em. MASK_Z2

     call dm_destroy(tmpw, ierr)

  else if (idx == 6) then

     ! dx_3d = repmat(dx, 1, 1, kb); 
     ! tmpdx1 = dx_3d(im, :, :);
     ! tmpdx2 = dx_3d(imm1, :, :);

     ! 
     ! tmpu1x = (2.e0 * u(im,:,:) * dti) ./ (tmpdx1 + tmpdx2);
     ! A1x = tmpu1x;
     ! A2x = tmpu1x;
     ! 
     ! A1x(tmpu1x > 0.e0) = 0.e0;
     ! A1x(tmpu1x <= 0.e0) = 1.e0;
     ! A2x(tmpu1x > 0.e0) = 1.e0;
     ! A2x(tmpu1x <= 0.e0) = 0.e0;
     ! 
     ! tmpq21x = q2(im, :, :) - (tmpu1x .* (small - q2(im, :, :)));
     ! tmpq2l1x = q2l(im, :, :) - (tmpu1x .* (small - q2l(im, :, :)));
     ! 
     ! tmpq22x = q2(im, :, :) - (tmpu1x .* (q2(im, :, :) - q2(imm1, :, :)));
     ! tmpq2l2x = q2l(im, :, :) - (tmpu1x .* (q2l(im, :, :) - q2l(imm1, :, :)));         
     ! 
     ! uf(im, :, :) = A1x(:,:,:) .* tmpq21x(:,:,:)+ A2x(:,:,:) .* tmpq22x(:,:,:);
     ! vf(im, :, :) = A1x(:,:,:) .* tmpq2l1x(:,:,:)+ A2x(:,:,:) .* tmpq2l2x(:,:,:);
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     ! 
     ! tmpdx1 = dx_3d(1, :, :);
     ! tmpdx2 = dx_3d(2, :, :);
     ! 
     ! tmpu1x = (2.e0 .* u(2, :, :) .* dti) ./ (tmpdx1 + tmpdx2);
     ! A1x = tmpu1x;
     ! A2x = tmpu1x;
     ! 
     ! A1x(tmpu1x >= 0.e0) = 1.e0;
     ! A1x(tmpu1x < 0.e0) = 0.e0;
     ! A2x(tmpu1x >= 0.e0) = 0.e0;
     ! A2x(tmpu1x < 0.e0) = 1.e0;
     ! 
     ! tmpq21x = q2(1, :, :) - tmpu1x .* (q2(1, :, :) - small);
     ! tmpq2l1x = q2l(1, :, :) - tmpu1x .* (q2l(1, :, :) - small);
     ! 
     ! tmpq22x = q2(1, :, :) - tmpu1x .* (q2(2, :, :) - q2(1, :, :));
     ! tmpq2l2x = q2l(1, :, :) - tmpu1x .* (q2l(2, :, :) - q2l(1, :, :));
     ! 
     ! uf(1, :, :) = A1x(1,:,:) .* tmpq21x(1,:,:) + A2x(1,:,:) .* tmpq22x(1,:,:);
     ! vf(1, :, :) = A1x(1,:,:) .* tmpq2l1x(1,:,:) + A2x(1,:,:) .* tmpq2l2x(1,:,:);
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! dy_3d = repmat(dy, 1, 1, kb);
     ! tmpdy1 = dy_3d(:, jm, :);
     ! tmpdy2 = dy_3d(:, jmm1, :);
     ! 
     ! tmpu1y = (2.e0 .* v(:, jm, :) .* dti) ./ (tmpdy1 + tmpdy2);
     ! A1y = tmpu1y;
     ! A2y = tmpu1y;
     ! 
     ! A1y(tmpu1y > 0.e0) = 0.e0;
     ! A1y(tmpu1y <= 0.e0) = 1.e0;
     ! A2y(tmpu1y > 0.e0) = 1.e0;
     ! A2y(tmpu1y <= 0.e0) = 0.e0;
     ! 
     ! tmpq21y = q2(:, jm, :) - tmpu1y .* (small - q2(:, jm, :));
     ! tmpq2l1y = q2l(:, jm, :) - tmpu1y .* (small - q2l(:, jm, :));
     ! 
     ! tmpq22y = q2(:, jm, :) - tmpu1y .* (q2(:, jm, :) - q2(:, jmm1, :));
     ! tmpq2l2y = q2l(:, jm, :) - tmpu1y .* (q2l(:, jm, :) - q2l(:, jmm1, :));
     ! 
     ! uf(:, jm, :) = tmpq21y(:,1,:) .* A1y(:,1,:) + tmpq22y(:,1,:) .* A2y(:,1,:);
     ! vf(:, jm, :) = tmpq2l1y(:,1,:) .* A1y(:,1,:) + tmpq2l2y(:,1,:) .* A2y(:,1,:);
     !   
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! tmpdy1 = dy_3d(:, 1, :);
     ! tmpdy2 = dy_3d(:, 2, :);
     ! 
     ! tmpu1y = (2.e0 .* v(:, 2, :) .* dti) ./ (tmpdy1 + tmpdy2);
     ! A1y = tmpu1y;
     ! A2y = tmpu1y;
     ! 
     ! A1y(tmpu1y >= 0.e0) = 1.e0;
     ! A1y(tmpu1y < 0.e0) = 0.e0;
     ! A2y(tmpu1y >= 0.e0) = 0.e0;
     ! A2y(tmpu1y < 0.e0) = 1.e0;
     ! 
     ! tmpq21y = q2(:, 1, :) - tmpu1y .* (q2(:, 1, :) - small);
     ! tmpq2l1y = q2l(:, 1, :) - tmpu1y .* (q2l(:, 1, :) - small);
     ! 
     ! tmpq22y = q2(:, 1, :) - tmpu1y .* (q2(:, 2, :) - q2(:, 1, :));
     ! tmpq2l2y = q2l(:, 1, :) - tmpu1y .* (q2l(:, 2, :) - q2l(:, 1, :));
     ! 
     ! uf(:, 1, :) = tmpq21y(:,1,:) .* A1y(:,1,:) + tmpq22y(:,1,:) .* A2y(:,1,:);
     ! vf(:, 1, :) = tmpq2l1y(:,1,:) .* A1y(:,1,:) + tmpq2l2y(:,1,:) .* A2y(:,1,:);
     ! 
     ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 
     ! ssmall = 1.e-10;
     ! tmpfsm = repmat(fsm, 1, 1, kb);
     ! 
     ! uf = uf .* tmpfsm + ssmall;
     ! vf = vf .* tmpfsm + ssmall;            

     dx_3d = dm_rep(dx, 1, 1, kb)

     tmp = (2.e0 * dti * (u .em. MASK_X2)) .ed. ((dx_3d .em. MASK_X2) + (shift(dx_3d, 1, -1) .em. MASK_X2))

     filter1 = (tmp <= 0.e0)
     filter2 = (tmp > 0.e0)

     filter1 = filter1 .em. MASK_X2
     filter2 = filter2 .em. MASK_X2

     tmp1 = q2 .em. MASK_X2 - (tmp .em. ((small - q2) .em. MASK_X2))
     tmp2 = q2l .em. MASK_X2 - (tmp .em. ((small - q2l) .em. MASK_X2))

     tmp3 = q2 .em. MASK_X2 - (tmp .em. (q2 .em. MASK_X2 - shift(q2, 1, -1) .em. MASK_X2))
     tmp4 = q2l .em. MASK_X2 - (tmp .em. (q2l .em. MASK_X2 - shift(q2l, 1, -1) .em. MASK_X2))

     uf = uf .em. REV_MASK_X2 + tmp1 .em. filter1 + tmp3 .em. filter2
     vf = vf .em. REV_MASK_X2 + tmp2 .em. filter1 + tmp4 .em. filter2


     tmp = (2.e0 * dti * (shift(u, 1, 1) .em. MASK_X1)) .ed. ((dx_3d .em. MASK_X1) + (shift(dx_3d, 1, 1) .em. MASK_X1))

     filter1 = (tmp >= 0.e0)
     filter2 = (tmp < 0.e0)

     filter1 = filter1 .em. MASK_X1
     filter2 = filter2 .em. MASK_X1

     tmp1 = q2 .em. MASK_X1 - tmp .em. ((q2 - small) .em. MASK_X1)
     tmp2 = q2l .em. MASK_X1 - tmp .em. ((q2l -small) .em. MASK_X1)

     tmp3 = q2 .em. MASK_X1 - tmp .em. (shift(q2, 1, 1) .em. MASK_X1 - q2 .em. MASK_X1)
     tmp4 = q2l .em. MASK_X1 - tmp .em. (shift(q2l, 1, 1) .em. MASK_X1 - q2l .em. MASK_X1)

     uf = uf .em. REV_MASK_X1 + tmp1 .em. filter1 + tmp3 .em. filter2
     vf = vf .em. REV_MASK_X1 + tmp2 .em. filter1 + tmp4 .em. filter2


     dy_3d = dm_rep(dy, 1, 1, kb)

     tmp = (2.e0 * dti * (v .em. MASK_Y2)) .ed. ((dy_3d .em. MASK_Y2) + &
          (shift(dy_3d, 2, -1) .em. MASK_Y2))

     filter1 = (tmp <= 0.e0)
     filter2 = (tmp > 0.e0)

     filter1 = filter1 .em. MASK_Y2
     filter2 = filter2 .em. MASK_Y2

     tmp1 = q2 .em. MASK_Y2 - (tmp .em. ((small - q2) .em. MASK_Y2))
     tmp2 = q2l .em. MASK_Y2 - (tmp .em. ((small - q2l) .em. MASK_Y2))

     tmp3 = q2 .em. MASK_Y2 - (tmp .em. (q2 .em. MASK_Y2 - shift(q2, 2, -1) .em. MASK_Y2))
     tmp4 = q2l .em. MASK_Y2 - (tmp .em. (q2l .em. MASK_Y2 - shift(q2l, 2, -1) .em. MASK_Y2))

     uf = uf .em. REV_MASK_Y2 + tmp1 .em. filter1 + tmp3 .em. filter2
     vf = vf .em. REV_MASK_Y2 + tmp2 .em. filter1 + tmp4 .em. filter2


     tmp = (2.e0 * dti * (shift(v, 2, 1) .em. MASK_Y1)) .ed. &
          ((dy_3d .em. MASK_Y1) + (shift(dy_3d, 2, 1) .em. MASK_Y1))

     filter1 = (tmp >= 0.e0)
     filter2 = (tmp < 0.e0)

     filter1 = filter1 .em. MASK_Y1
     filter2 = filter2 .em. MASK_Y1

     tmp1 = q2 .em. MASK_Y1 - tmp  .em. ((q2 - small) .em. MASK_Y1)
     tmp2 = q2l .em. MASK_Y1 - tmp .em. ((q2l -small) .em. MASK_Y1)

     tmp3 = q2 .em. MASK_Y1 - tmp .em. (shift(q2, 2, 1) .em. MASK_Y1 - q2 .em. MASK_Y1)
     tmp4 = q2l .em. MASK_Y1 - tmp .em. (shift(q2l, 2, 1) .em. MASK_Y1 - q2l .em. MASK_Y1)

     uf = uf .em. REV_MASK_Y1 + tmp1 .em. filter1 + tmp3 .em. filter2
     vf = vf .em. REV_MASK_Y1 + tmp2 .em. filter1 + tmp4 .em. filter2

     uf = 1.e-10 + uf .em. dm_rep(fsm, 1, 1, kb)
     vf = 1.e-10 + vf .em. dm_rep(fsm, 1, 1, kb)

     call dm_destroy(dx_3d, ierr)
     call dm_destroy(dy_3d, ierr)
     call dm_destroy(tmp, ierr)
     call dm_destroy(tmp1, ierr)
     call dm_destroy(tmp2, ierr)
     call dm_destroy(tmp3, ierr)
     call dm_destroy(tmp4, ierr)
     call dm_destroy(filter1, ierr)
     call dm_destroy(filter2, ierr)
  end if

  deallocate(arr_im, arr_jm, seq_im, seq_jm, arr_imm1, arr_jmm1, arr_imm2, arr_jmm2)
  
end subroutine new_bcond
