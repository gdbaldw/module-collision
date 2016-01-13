# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/Makefile.template,v 1.3 2007/02/26 00:34:32 masarati Exp $
###############################################################################
#
# Add module dependencies here; e.g., if in dir 'module-X' there's 
# a file 'myfunc.cc' write:
#
# MODULE_DEPENDENCIES=myfunc.o
#
# Then put this file in directory 'module-X' as 'Makefile.inc'
#
###############################################################################

MODULE_DEPENDENCIES=intersect.lo
MODULE_LINK=-lfcl
