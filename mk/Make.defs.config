# -*- Mode: Makefile;  -*-

## Set the configuration string used in filenames and directories
## This file should not need to be edited by users who are trying to compile
## Developers should change it only when they are adding 
## a new module or feature that users can toggle.

## Usage:
##  When this file is `include'd in another makefile, the following
##  variables must already be defined:
##     ENV == name of the "ENV" variable that is a standin for compiler options
##     USE_WWM == include the wind-wave module 
##     DEBUG == TRUE for symbol table, FALSE for no symbol table
##     OPT == FALSE no optimization
##            TRUE optimization, asserts, and initial setVal
##            HIGH optimization, no asserts, and initialize to zero
##     ENV == name of the user environment


makefiles+=Make.defs.config

# If OPT isn't set, set it to the opposite of DEBUG
ifeq ($(OPT),)
  ifeq ($(DEBUG),yes)
    OPT := no
  endif
  ifeq ($(DEBUG),no)
    OPT := yes
  endif
endif

# Set USE_SETVAL to TRUE unless OPT=HIGH
ifeq ($(USE_SETVAL),)
  USE_SETVAL := TRUE
  ifeq ($(OPT),HIGH)
    USE_SETVAL := no
  endif
endif



# these variables specify pieces of the configuration string
_vers := svn

ifneq ($(ENV),)
  _env  := .$(ENV)
endif

_debug  := $(subst no,,$(subst yes,.DEBUG,$(USE_DEBUG)))
_opt    := $(subst no,,$(subst yes,.OPT,$(subst HIGH,.OPTHIGH,$(OPT))))
_wwm    := $(subst no,,$(subst yes,.WWM,$(USE_WWM)))
_petsc  := $(subst no,,$(subst yes,.PETSC,$(USE_PETSC)))
_sed    := $(subst no,,$(subst yes,.SED,$(USE_SED)))
_sed2d  := $(subst no,,$(subst yes,.SED2D,$(USE_SED2D)))
_eco    := $(subst no,,$(subst yes,.ECO,$(USE_ECO)))
_icm    := $(subst no,,$(subst yes,.ICM,$(USE_ICM)))
_napzd  := $(subst no,,$(subst yes,.NAPZD,$(USE_NAPZD)))
_timer  := $(subst no,,$(subst yes,.TIMER,$(USE_TIMER)))


# create a base for all of the other config strings
config  := $(shell echo $(_vers)$(_env)$(_wwm)$(_petsc)$(_sed2d)$(_sed)$(_eco)$(_icm)$(_napzd)$(_timer)$(_debug)$(_opt) | sed -e 's/ //g' -e 's/	//g')
