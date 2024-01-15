# Hey Emacs, this is a -*- makefile -*-
#----------------------------------------------------------------------------
#
# Makefile for ChipWhisperer SimpleSerial-Masked-Dilithium Program
#
#----------------------------------------------------------------------------
# On command line:
#
# make all = Make software.

# make clean = Clean out built project files.
#
# make coff = Convert ELF to AVR COFF.
#
# make extcoff = Convert ELF to AVR Extended COFF.
#
# make program = Download the hex file to the device, using avrdude.
#                Please customize the avrdude settings below first!
#
# make debug = Start either simulavr or avarice as specified for debugging,
#              with avr-gdb or avr-insight as the front end for debugging.
#
# make filename.s = Just compile filename.c into the assembler code only.
#
# make filename.i = Create a preprocessed source file for use in submitting
#                   bug reports to the GCC project.
#
# To rebuild project do "make clean" then "make all".
#----------------------------------------------------------------------------

# Target file name (without extension).
# This is the base name of the compiled .hex file.
TARGET = kyber-masked-cmp
VID = 0x2b3e
PID ?= 0xace2
CW_OPENOCD_CFG ?= ~/chipwhisperer/openocd/cw_openocd.cfg
# NOTE: set TARGET_CFG to target/stm32{f3x,f4x,l5x}.cfg if needed
TARGET_CFG ?= target/stm32f4x.cfg
SS_VER ?= SS_VER_2_1
# optimization level
OPT ?= 3

SRCDIR = src
CDEFS += -DCW -DCOMPARE_USE_A2A
CFLAGS += -I${SRCDIR} -Wall -Wextra
CFLAGS += -g
CFLAGS += -DCW_MEASURE

# CFLAGS += {-DNSHARES>=2}
CFLAGS += -DNSHARES=4
# CFLAGS += {-DKYBER, -DSABER}
CFLAGS += -DKYBER
# CFLAGS += {-DGF, -DARITH, -DSIMPLE, -DHYBRIDSIMPLE}
CFLAGS += -DGF
# CFLAGS += {-NTESTS>=1}
#CFLAGS += -DNTESTS=10000
CFLAGS += -DNTESTS=3

# List C source files here.
# Header files (.h) are automatically pulled in.
#CPPSRC += simpleserial-aes.cpp
SRC += kyber-masked-cmp.c
SRC += $(wildcard ${SRCDIR}/*.c)

# -----------------------------------------------------------------------------
EXTRA_OPTS = NO_EXTRA_OPTS
CFLAGS += -D$(EXTRA_OPTS)


ifeq ($(CRYPTO_TARGET),)
  ${info No CRYPTO_TARGET passed - defaulting to NONE}
  CRYPTO_TARGET = NONE
endif

# NOTE: set to CW308_STM32{F3,F4,L5} if needed
ifeq ($(PLATFORM),)
  ${info No PLATFORM passed - defaulting to CW308_STM32F4}
  PLATFORM = CW308_STM32F4
endif

${info Building for platform ${PLATFORM} with CRYPTO_TARGET=$(CRYPTO_TARGET)}

#Add simpleserial project to build
include ../simpleserial/Makefile.simpleserial

FIRMWAREPATH = ../.
include $(FIRMWAREPATH)/Makefile.inc

# quick and dirty
$(shell mkdir -p objdir-${PLATFORM}/${SRCDIR})
