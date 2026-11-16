#!/bin/bash
# Marge: single node, no InfiniBand card to choose; the GPU is set in main.f90
LRANK=${OMPI_COMM_WORLD_LOCAL_RANK:-$SLURM_LOCALID}
GRANK=${OMPI_COMM_WORLD_RANK:-$SLURM_PROCID}
echo "Rank $GRANK local $LRANK on $(hostname)"
# "$@" keeps each argument intact
exec "$@"
