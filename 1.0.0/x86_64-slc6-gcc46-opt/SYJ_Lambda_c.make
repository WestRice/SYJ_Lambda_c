#-- start of make_header -----------------

#====================================
#  Library SYJ_Lambda_c
#
#   Generated Mon Aug 20 12:59:12 2018  by suyj
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_SYJ_Lambda_c_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_SYJ_Lambda_c_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_SYJ_Lambda_c

SYJ_Lambda_c_tag = $(tag)

#cmt_local_tagfile_SYJ_Lambda_c = $(SYJ_Lambda_c_tag)_SYJ_Lambda_c.make
cmt_local_tagfile_SYJ_Lambda_c = $(bin)$(SYJ_Lambda_c_tag)_SYJ_Lambda_c.make

else

tags      = $(tag),$(CMTEXTRATAGS)

SYJ_Lambda_c_tag = $(tag)

#cmt_local_tagfile_SYJ_Lambda_c = $(SYJ_Lambda_c_tag).make
cmt_local_tagfile_SYJ_Lambda_c = $(bin)$(SYJ_Lambda_c_tag).make

endif

include $(cmt_local_tagfile_SYJ_Lambda_c)
#-include $(cmt_local_tagfile_SYJ_Lambda_c)

ifdef cmt_SYJ_Lambda_c_has_target_tag

cmt_final_setup_SYJ_Lambda_c = $(bin)setup_SYJ_Lambda_c.make
cmt_dependencies_in_SYJ_Lambda_c = $(bin)dependencies_SYJ_Lambda_c.in
#cmt_final_setup_SYJ_Lambda_c = $(bin)SYJ_Lambda_c_SYJ_Lambda_csetup.make
cmt_local_SYJ_Lambda_c_makefile = $(bin)SYJ_Lambda_c.make

else

cmt_final_setup_SYJ_Lambda_c = $(bin)setup.make
cmt_dependencies_in_SYJ_Lambda_c = $(bin)dependencies.in
#cmt_final_setup_SYJ_Lambda_c = $(bin)SYJ_Lambda_csetup.make
cmt_local_SYJ_Lambda_c_makefile = $(bin)SYJ_Lambda_c.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)SYJ_Lambda_csetup.make

#SYJ_Lambda_c :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'SYJ_Lambda_c'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = SYJ_Lambda_c/
#SYJ_Lambda_c::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

SYJ_Lambda_clibname   = $(bin)$(library_prefix)SYJ_Lambda_c$(library_suffix)
SYJ_Lambda_clib       = $(SYJ_Lambda_clibname).a
SYJ_Lambda_cstamp     = $(bin)SYJ_Lambda_c.stamp
SYJ_Lambda_cshstamp   = $(bin)SYJ_Lambda_c.shstamp

SYJ_Lambda_c :: dirs  SYJ_Lambda_cLIB
	$(echo) "SYJ_Lambda_c ok"

#-- end of libary_header ----------------

SYJ_Lambda_cLIB :: $(SYJ_Lambda_clib) $(SYJ_Lambda_cshstamp)
	@/bin/echo "------> SYJ_Lambda_c : library ok"

$(SYJ_Lambda_clib) :: $(bin)Sigma0PionPi0.o $(bin)Sigma0PionEta.o $(bin)Entries.o $(bin)Load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(SYJ_Lambda_clib) $?
	$(lib_silent) $(ranlib) $(SYJ_Lambda_clib)
	$(lib_silent) cat /dev/null >$(SYJ_Lambda_cstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(SYJ_Lambda_clibname).$(shlibsuffix) :: $(SYJ_Lambda_clib) $(SYJ_Lambda_cstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" SYJ_Lambda_c $(SYJ_Lambda_c_shlibflags)

$(SYJ_Lambda_cshstamp) :: $(SYJ_Lambda_clibname).$(shlibsuffix)
	@if test -f $(SYJ_Lambda_clibname).$(shlibsuffix) ; then cat /dev/null >$(SYJ_Lambda_cshstamp) ; fi

SYJ_Lambda_cclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)Sigma0PionPi0.o $(bin)Sigma0PionEta.o $(bin)Entries.o $(bin)Load.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
SYJ_Lambda_cinstallname = $(library_prefix)SYJ_Lambda_c$(library_suffix).$(shlibsuffix)

SYJ_Lambda_c :: SYJ_Lambda_cinstall

install :: SYJ_Lambda_cinstall

SYJ_Lambda_cinstall :: $(install_dir)/$(SYJ_Lambda_cinstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(SYJ_Lambda_cinstallname) :: $(bin)$(SYJ_Lambda_cinstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(SYJ_Lambda_cinstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(SYJ_Lambda_cinstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(SYJ_Lambda_cinstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(SYJ_Lambda_cinstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(SYJ_Lambda_cinstallname) $(install_dir)/$(SYJ_Lambda_cinstallname); \
	      echo `pwd`/$(SYJ_Lambda_cinstallname) >$(install_dir)/$(SYJ_Lambda_cinstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(SYJ_Lambda_cinstallname), no installation directory specified"; \
	  fi; \
	fi

SYJ_Lambda_cclean :: SYJ_Lambda_cuninstall

uninstall :: SYJ_Lambda_cuninstall

SYJ_Lambda_cuninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(SYJ_Lambda_cinstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(SYJ_Lambda_cinstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(SYJ_Lambda_cinstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(SYJ_Lambda_cinstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of dependencies ------------------
ifneq ($(MAKECMDGOALS),SYJ_Lambda_cclean)
ifneq ($(MAKECMDGOALS),uninstall)

$(bin)SYJ_Lambda_c_dependencies.make : $(use_requirements) $(cmt_final_setup_SYJ_Lambda_c)
	$(echo) "(SYJ_Lambda_c.make) Rebuilding $@"; \
	  $(build_dependencies) -out=$@ -start_all  -end_all $(includes) $(app_SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) -name=SYJ_Lambda_c $? -f=$(cmt_dependencies_in_SYJ_Lambda_c) -without_cmt

-include $(bin)SYJ_Lambda_c_dependencies.make

endif
endif

SYJ_Lambda_cclean ::
	$(cleanup_silent) \rm -rf $(bin)SYJ_Lambda_c_deps $(bin)SYJ_Lambda_c_dependencies.make
#-- end of dependencies -------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SYJ_Lambda_cclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Sigma0PionPi0.d

$(bin)$(binobj)Sigma0PionPi0.d :

$(bin)$(binobj)Sigma0PionPi0.o : $(cmt_final_setup_SYJ_Lambda_c)

$(bin)$(binobj)Sigma0PionPi0.o : $(src)Sigma0PionPi0.cxx
	$(cpp_echo) $(src)Sigma0PionPi0.cxx
	$(cpp_silent) $(cppcomp)  -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Sigma0PionPi0_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Sigma0PionPi0_cppflags) $(Sigma0PionPi0_cxx_cppflags)  $(src)Sigma0PionPi0.cxx
endif
endif

else
$(bin)SYJ_Lambda_c_dependencies.make : $(Sigma0PionPi0_cxx_dependencies)

$(bin)SYJ_Lambda_c_dependencies.make : $(src)Sigma0PionPi0.cxx

$(bin)$(binobj)Sigma0PionPi0.o : $(Sigma0PionPi0_cxx_dependencies)
	$(cpp_echo) $(src)Sigma0PionPi0.cxx
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Sigma0PionPi0_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Sigma0PionPi0_cppflags) $(Sigma0PionPi0_cxx_cppflags)  $(src)Sigma0PionPi0.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SYJ_Lambda_cclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Sigma0PionEta.d

$(bin)$(binobj)Sigma0PionEta.d :

$(bin)$(binobj)Sigma0PionEta.o : $(cmt_final_setup_SYJ_Lambda_c)

$(bin)$(binobj)Sigma0PionEta.o : $(src)Sigma0PionEta.cxx
	$(cpp_echo) $(src)Sigma0PionEta.cxx
	$(cpp_silent) $(cppcomp)  -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Sigma0PionEta_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Sigma0PionEta_cppflags) $(Sigma0PionEta_cxx_cppflags)  $(src)Sigma0PionEta.cxx
endif
endif

else
$(bin)SYJ_Lambda_c_dependencies.make : $(Sigma0PionEta_cxx_dependencies)

$(bin)SYJ_Lambda_c_dependencies.make : $(src)Sigma0PionEta.cxx

$(bin)$(binobj)Sigma0PionEta.o : $(Sigma0PionEta_cxx_dependencies)
	$(cpp_echo) $(src)Sigma0PionEta.cxx
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Sigma0PionEta_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Sigma0PionEta_cppflags) $(Sigma0PionEta_cxx_cppflags)  $(src)Sigma0PionEta.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SYJ_Lambda_cclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Entries.d

$(bin)$(binobj)Entries.d :

$(bin)$(binobj)Entries.o : $(cmt_final_setup_SYJ_Lambda_c)

$(bin)$(binobj)Entries.o : $(src)components/Entries.cxx
	$(cpp_echo) $(src)components/Entries.cxx
	$(cpp_silent) $(cppcomp)  -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Entries_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Entries_cppflags) $(Entries_cxx_cppflags) -I../src/components $(src)components/Entries.cxx
endif
endif

else
$(bin)SYJ_Lambda_c_dependencies.make : $(Entries_cxx_dependencies)

$(bin)SYJ_Lambda_c_dependencies.make : $(src)components/Entries.cxx

$(bin)$(binobj)Entries.o : $(Entries_cxx_dependencies)
	$(cpp_echo) $(src)components/Entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Entries_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Entries_cppflags) $(Entries_cxx_cppflags) -I../src/components $(src)components/Entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),SYJ_Lambda_cclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Load.d

$(bin)$(binobj)Load.d :

$(bin)$(binobj)Load.o : $(cmt_final_setup_SYJ_Lambda_c)

$(bin)$(binobj)Load.o : $(src)components/Load.cxx
	$(cpp_echo) $(src)components/Load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Load_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Load_cppflags) $(Load_cxx_cppflags) -I../src/components $(src)components/Load.cxx
endif
endif

else
$(bin)SYJ_Lambda_c_dependencies.make : $(Load_cxx_dependencies)

$(bin)SYJ_Lambda_c_dependencies.make : $(src)components/Load.cxx

$(bin)$(binobj)Load.o : $(Load_cxx_dependencies)
	$(cpp_echo) $(src)components/Load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(SYJ_Lambda_c_pp_cppflags) $(lib_SYJ_Lambda_c_pp_cppflags) $(Load_pp_cppflags) $(use_cppflags) $(SYJ_Lambda_c_cppflags) $(lib_SYJ_Lambda_c_cppflags) $(Load_cppflags) $(Load_cxx_cppflags) -I../src/components $(src)components/Load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: SYJ_Lambda_cclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(SYJ_Lambda_c.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

SYJ_Lambda_cclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library SYJ_Lambda_c
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)SYJ_Lambda_c$(library_suffix).a $(library_prefix)SYJ_Lambda_c$(library_suffix).s? SYJ_Lambda_c.stamp SYJ_Lambda_c.shstamp
#-- end of cleanup_library ---------------
