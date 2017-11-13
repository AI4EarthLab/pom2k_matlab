# pom2k_matlab
Matlab code for pom2k ocean model
通过比较pom08-wad 与pom2K差异，在Matlab code中添加干湿网格，
代码先由matlab_orgin复制到matlab_wad
比较方法：
pom08-wad是基于pom2k稍前的版本改造的，基础版本有少量差异要先进行一定修改才方便比较，
格式修改：pom2k是Fortran90 自由与固定兼容格式，pom08-wad是固定格式
首先把pom08-wad改为自由与固定兼容格式;
INCLUE修改：由于include文件名不同，用宏定义进行include
使用windiff进行文件比较，在不同的地方进行人工分析，用宏标记，标记为：
变量不同：ZWADMODIFY_V   real:: dummy##__line__=wadmodify
代码不同：ZWADMODIFY_C wadmodify=__line__
再用understand 进行分析wadmodify的引用情况
    Set [pom2008.f90,  261]        pom2008
    Set [pom2008.f90,  296]        pom2008
    Set [pom2008.f90,  387]        pom2008
    Set [pom2008.f90,  456]        pom2008
    Set [pom2008.f90,  461]        pom2008
    Set [pom2008.f90,  640]        pom2008
    Set [pom2008.f90,  671]        pom2008
    Set [pom2008.f90,  851]        pom2008
    Set [pom2008.f90,  862]        pom2008
    Set [pom2008.f90,  911]        pom2008
    Set [pom2008.f90,  947]        pom2008
    Set [pom2008.f90,  972]        pom2008
    Set [pom2008.f90,  988]        pom2008
    Set [pom2008.f90, 1007]        pom2008
    Set [pom2008.f90, 1152]        pom2008
    Set [pom2008.f90, 1220]        pom2008
    Set [pom2008.f90, 1312]        pom2008
    Set [pom2008.f90, 1351]        pom2008 OK    
    Set [pom2008.f90, 1461]        pom2008
    Set [pom2008.f90, 1524]        pom2008
    Set [pom2008.f90, 1668]        pom2008
    Set [pom2008.f90, 1703]        pom2008
    Set [pom2008.f90, 1748]        pom2008
    Set [pom2008.f90, 1867]        pom2008
    Set [pom2008.f90, 1879]        pom2008
    Set [pom2008.f90, 1960]        pom2008
    Set [pom2008.f90, 2975]        areas_masks
    Set [pom2008.f90, 3281]        bcond
    Set [pom2008.f90, 3293]        bcond
    Set [pom2008.f90, 3308]        bcond
    Set [pom2008.f90, 3318]        bcond
    Set [pom2008.f90, 3351]        bcond
    Set [pom2008.f90, 3362]        bcond
    Set [pom2008.f90, 3378]        bcond
    Set [pom2008.f90, 3388]        bcond
    Set [pom2008.f90, 4096]        box
    Set [pom2008.f90, 4137]        box
    Set [pom2008.f90, 4241]        dens
    Set [pom2008.f90, 4486]        file2ic
    Set [pom2008.f90, 4580]        file2ic
    Set [pom2008.f90, 4587]        file2ic
    Set [pom2008.f90, 4977]        profq
    Set [pom2008.f90, 4997]        profq
    Set [pom2008.f90, 5183]        profq
    Set [pom2008.f90, 6032]        seamount
    Set [pom2008.f90, 6204]        seamount
    Set [pom2008.f90, 6209]        seamount
    Set [pom2008.f90, 6274]        seamount
    Set [pom2008.f90, 6295]        seamount
根据分析进行修改
相关代码在matlab_wad/fcode