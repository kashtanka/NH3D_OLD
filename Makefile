#makefile for nh3d nonparallel
model_path = ~/NH3D_OLD/
source_path = ~/NH3D_OLD/source
obj_path = ~/NH3D_OLD/obj
exec = nh3d.out
fc = ifort
switch = -r8 -O2 #-check all -traceback 
obj = \
 $(obj_path)/modules.o   \
 $(obj_path)/nhad0626.o  \
 $(obj_path)/nhli0606.o  \
 $(obj_path)/eq_state.o  \
 $(obj_path)/moment.o    \
 $(obj_path)/ellipt.o    \
 $(obj_path)/tracers.o   \
 $(obj_path)/diffu_local.o     \
 $(obj_path)/dif_Noh.o   \
 $(obj_path)/closures.o  \
 $(obj_path)/bl_depth_sub.o \
 $(obj_path)/cliradlw.o  \
 $(obj_path)/cliradtfc.o \
 $(obj_path)/dst3_adv.o  \
 $(obj_path)/fct2d.o     \
 $(obj_path)/hsmoot2.o   \
 $(obj_path)/hsmoot_6_lim.o \
 $(obj_path)/hsmop2.o    \
 $(obj_path)/hsmop_6_lim.o  \
 $(obj_path)/impur.o     \
 $(obj_path)/inrefs_sub.o \
 $(obj_path)/LS_closure.o \
 $(obj_path)/microphysics.o \
 $(obj_path)/out6_sub.o   \
 $(obj_path)/perturb_profile.o \
 $(obj_path)/surface_layer.o \
 $(obj_path)/wind_perturb.o \
#
$(exec) : $(obj)
	$(fc) -o $(exec) $(switch) $(obj)
#MODULES:
$(obj_path)/modules.o : $(source_path)/modules.f
	cd $(obj_path)/ && $(fc) -c $(switch) $(source_path)/modules.f && cd $(model_path)
#OTHER SOURCE
$(obj_path)/%.o : $(source_path)/%.for
	cd $(obj_path)/ && $(fc) -c $(switch) $< && cd $(model_path) 
$(obj_path)/%.o : $(source_path)/%.f90
	cd $(obj_path)/ && $(fc) -c $(switch) $< && cd $(model_path)
$(obj_path)/%.o : $(source_path)/%.f
	cd $(obj_path)/ && $(fc) -c $(switch) $< && cd $(model_path)
clean :
	rm $(obj) $(exec)
