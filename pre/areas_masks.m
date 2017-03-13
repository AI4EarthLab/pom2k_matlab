function areas_masks()
global im imm1 jm jmm1 kb dx dy h art aru arv fsm dum dvm art_3d aru_3d arv_3d fsm_3d dum_3d dvm_3d;

% Calculate areas of "t" and "s" cells:
art = dx .* dy;
aru(2:im,2:jm)=.25e0.*(dx(2:im,2:jm)+dx(1:imm1,2:jm)).*(dy(2:im,2:jm)+dy(1:imm1,2:jm));
arv(2:im,2:jm)=.25e0.*(dx(2:im,2:jm)+dx(2:im,1:jmm1)).*(dy(2:im,2:jm)+dy(2:im,1:jmm1));
aru(1,:)=aru(2,:);      aru(:,1)=aru(:,2);
arv(1,:)=arv(2,:);      arv(:,1)=arv(:,2);
art_3d=repmat(art,1,1,kb);  aru_3d=repmat(aru,1,1,kb);  arv_3d=repmat(arv,1,1,kb);

% Initialise and set up free surface mask:
fsm(h>1.0)=1.0;

dum(2:im,2:jm)=fsm(2:im,2:jm) .* fsm(1:im-1,2:jm);
dvm(2:im,2:jm)=fsm(2:im,2:jm) .* fsm(2:im,1:jm-1);
dum_3d=repmat(dum,1,1,kb);  dvm_3d=repmat(dvm,1,1,kb);  fsm_3d=repmat(fsm,1,1,kb);