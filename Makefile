TOPDIR = .
SUBDIRS := $(shell find * -maxdepth 0 -type d | grep -v SpNMerge)

#Exclude these directories
EXCLUDEDIR := Analysis Doc Fortran HiRep.xcodeproj TestProgramUtils
SUBDIRS := $(filter-out $(EXCLUDEDIR), $(SUBDIRS))

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules


