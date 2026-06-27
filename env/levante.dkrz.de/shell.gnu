# Levante (DKRZ) — GNU toolchain. Portability build (not the bit-identity anchor).
export LC_ALL=en_US.UTF-8
export CPU_MODEL=AMD_EPYC_ZEN3

module --force purge

module load gcc/11.2.0-gcc-11.2.0
module load openmpi/4.1.2-gcc-11.2.0
export FC=mpif90 CC=mpicc CXX=mpicxx

module load netcdf-c/4.8.1-gcc-11.2.0
module load netcdf-fortran/4.5.3-gcc-11.2.0
module load git

ulimit -s unlimited 2>/dev/null || true   # best-effort (see shell.intel): on a SLURM compute node
ulimit -c 0 2>/dev/null || true           # the hard stack may be capped; raising it errors under set -e

export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl=self
export OMPI_MCA_osc="pt2pt"
export UCX_IB_ADDR_TYPE=ib_global
export OMPI_MCA_coll="^ml,hcoll"
export OMPI_MCA_coll_hcoll_enable="0"
export HCOLL_ENABLE_MCAST_ALL="0"
export HCOLL_MAIN_IB=mlx5_0:1
export UCX_NET_DEVICES=mlx5_0:1
export UCX_TLS=mm,knem,cma,dc_mlx5,dc_x,self
export UCX_UNIFIED_MODE=y
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_io="romio321"
export UCX_HANDLE_ERRORS=bt
