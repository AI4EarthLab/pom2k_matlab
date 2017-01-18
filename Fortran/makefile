
PNETCDF=${ESM_SOFT}/pnetcdf/

CFLAGS	         =  
FFLAGS	         =-Wno-tabs -I ${PNETCDF}/include  -J ${DM} -g
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = pom2k *.o *.mod *.nc
NP               = 1
OBJ		 = input.o grid.o  file2ic.o proft.o profu.o  \
		   vertvl.o advave.o baropg.o advct.o advt1.o \
                   advu.o advv.o dense.o depth.o advq.o profq.o profv.o bcond.o

OBJMAIN		 = ${OBJ} pom2k.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  -${FLINKER} -o pom2k ${OBJMAIN} ${PETSC_KSP_LIB} \
      ${PNETCDF}/lib/libpnetcdf.a -J ${DM} ${DM}/libdm.a

run : main
	-@./pom2k

# small:
# 	make clean
# 	make main 
# 	-@${MPIEXEC} -n 4 ./main -m 3 -n 2 -k 3 -ep 3.1 -debug -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type jacobi



#include ${PETSC_DIR}/lib/petsc/conf/test
