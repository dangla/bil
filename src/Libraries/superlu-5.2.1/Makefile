#
#
SHELL = /bin/sh


dFILES   = ${wildcard ${addsuffix d*.c,./SRC/}}


ALLAUX  = superlu_timer.c util.c memory.c get_perm_c.c mmd.c \
          sp_coletree.c sp_preorder.c sp_ienv.c relax_snode.c \
          heap_relax_snode.c colamd.c \
          ilu_relax_snode.c ilu_heap_relax_snode.c mark_relax.c \
          mc64ad.c qselect.c input_error.c smach.c dmach.c

DLUSRC = \
        dgssv.c dgssvx.c \
        dsp_blas2.c dsp_blas3.c dgscon.c \
        dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c  \
        dgsrfs.c dgstrf.c dgstrs.c dcopy_to_ucol.c \
        dsnode_dfs.c dsnode_bmod.c dpanel_dfs.c dpanel_bmod.c \
        dreadhb.c dreadrb.c dreadtriple.c \
        dcolumn_dfs.c dcolumn_bmod.c dpivotL.c dpruneL.c \
        dmemory.c dutil.c dmyblas2.c \
        dgsisx.c dgsitrf.c dldperm.c \
        ilu_ddrop_row.c ilu_dsnode_dfs.c \
        ilu_dcolumn_dfs.c ilu_dpanel_dfs.c ilu_dcopy_to_ucol.c \
        ilu_dpivotL.c ddiagonal.c dlacon2.c
        ## dgstrsL.c dgstrsU.c

FILES = ${DLUSRC} ${ALLAUX}

test:
	@echo "#if defined(__cplusplus)"
	@echo "  extern "C" {"
	@echo "#endif"
	@echo
	@for i in ${FILES}; do \
	  echo "#include \"./SRC/$$i\"" ; \
	done
	@echo
	@echo "#if defined(__cplusplus)"
	@echo "  }"
	@echo "#endif"
