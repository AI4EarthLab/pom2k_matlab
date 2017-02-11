function [tb,sb,tclim,sclim,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs]=init_ts(grid,h,zz)
    global im jm kb kbm1
    global tbias sbias

    M0=zeros(im,jm,kb);
    cor     = Field(M0, grid, 3);
    tb      = Field(M0, grid, 3);
    sb      = Field(M0, grid, 3);
    tclim   = Field(M0, grid, 3);
    sclim   = Field(M0, grid, 3);
    
    tbe     = Field(M0, grid, 3);
    tbw     = Field(M0, grid, 3);
    sbe     = Field(M0, grid, 3);
    sbw     = Field(M0, grid, 3);
    tbn     = Field(M0, grid, 3);
    tbs     = Field(M0, grid, 3);
    sbn     = Field(M0, grid, 3);
    sbs     = Field(M0, grid, 3);
  
    % Set initial conditions:
    for k=1:kbm1
        tb(:,:,k)=5.0+15.0*exp(zz(k)*h(:,:)/1000.0)-tbias;
        sb(:,:,k)=35.0-sbias;
        tclim(:,:,k)=tb(:,:,k);
        sclim(:,:,k)=sb(:,:,k);
    end 
    tbe(:,1:kbm1) = tb(im,:,1:kbm1);
    tbw(:,1:kbm1) = tb(1,:,1:kbm1);
    sbe(:,1:kbm1) = sb(im,:,1:kbm1);
    sbw(:,1:kbm1) = sb(1,:,1:kbm1);
    tbn(:,1:kbm1) = tb(:,jm,1:kbm1);
    tbs(:,1:kbm1) = tb(:,1,1:kbm1);
    sbn(:,1:kbm1) = sb(:,jm,1:kbm1);
    sbs(:,1:kbm1) = sb(:,1,1:kbm1);
end