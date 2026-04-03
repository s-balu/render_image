# =============================================================================
# imbas_renderer — Config.mk
# =============================================================================

ENABLE_OPENMP ?= 1	# OpenMP on? 1 yes, 0 no
ENABLE_MPI    ?= 0	# MPI on? 1 yes, 0 no

OPT += -DCIC		# Include CIC
OPT += -DCIC_STRIP_HEIGHT=64
OPT += -DMAX_H_DEPOSIT_PX=8
OPT += -DMIN_RHO_FRAC=0.1
OPT += -DIMAGE_DIMENSIONX=1024
OPT += -DIMAGE_DIMENSIONY=1024
#OPT += -DMAX_H_VOXELS=8
#OPT += -DLONG_IDS
#OPT += -DVELOCITIES
#OPT += -DIMAGE
#OPT += -DDEBUG
#OPT += -DSCATTER_DECOMPOSITION
#OPT += -DNONPERIODIC
#OPT += -DKERNEL_SMOOTHING
#OPT += -DSMOOTH_GRID=512

# ---------------------------------------------------------------------------
# Target system  (uncomment exactly one)
# ---------------------------------------------------------------------------
SYSTEM ?= MacBook
#SYSTEM ?= MacPro
#SYSTEM ?= OzSTAR
#SYSTEM ?= Setonix
#SYSTEM ?= Naranjo

