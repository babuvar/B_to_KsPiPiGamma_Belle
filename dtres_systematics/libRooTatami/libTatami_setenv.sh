############################################
# sample .bashrc by S.Nishida (2016/08/22) #
############################################
#
# *** Please copy following lines to your ~/.bashrc to set up
# *** Belle analysis environment
#
# Note: USE_GRAND_REPROCESS_DATA is now forcely set to 1 at 
#       /sw/belle/local/etc/bashrc

#####################
# Belle Environment #
#####################
export BELLE_LEVEL=b20090127_0910
export BELLE_DEBUG=opt
export BELLE_MSG_MAX_SHLVL=1
# use i386 version rather than RHEL4 64bit version for exp5x MC prod
export USE_I386_EXP5X=
# avoid glibc detected error
export MALLOC_CHECK_=0

#
# *** by default, ROOT version 5.34.36 is used. If you want to use
# *** different version, please uncomment and modify below
#
# export ROOTSYS=/sw/belle/cern/root_v5.xx.xx

if [ -f /sw/belle/local/etc/bashrc_general ]; then
    . /sw/belle/local/etc/bashrc_general
fi

# setup ROOT version
export ROOTSYS=/home/belle/varghese/IPHC/old-root/root-tatami
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:.:$LD_LIBRARY_PATH
export MANPATH=$MANPATH:$ROOTSYS/man

# setup BELLE_FLC
#export BELLE_FLC=/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way
export BELLE_FLC=`pwd`
export LD_LIBRARY_PATH=$BELLE_FLC/lib:$LD_LIBRARY_PATH



