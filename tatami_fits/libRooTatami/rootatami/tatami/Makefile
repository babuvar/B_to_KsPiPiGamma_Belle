MAKE = gmake

srcdir := $(subst config/../,,.)
top_srcdir := $(subst config/../,,..)
VPATH := $(subst config/../,,.)

BELLE_CONFIG_DIR = /sw/belle/belle/b20090127_0910/src/config

PACKAGE = tatami
VERSION = 1.0

packagedir = .

PACKAGE_HEADER_DIR = $(PACKAGE)
export PACKAGE_HEADER_DIR
BELLE_ROUTE = src/anal/tatami/
export BELLE_ROUTE

SUBDIRS = src tatami
MAINDIRS = mains 

include $(BELLE_CONFIG_DIR)/Makefile.subdirs

