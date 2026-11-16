#!/bin/bash
# Vega GPU nodes: InfiniBand card of each local rank, as in MHIT36 (mlx5_0, mlx5_1);
# ranks 2,3 assumed on the second card (check with nvidia-smi topo -m)
# local rank: from mpirun (OMPI_*) or from srun (SLURM_*)
LRANK=${OMPI_COMM_WORLD_LOCAL_RANK:-$SLURM_LOCALID}
GRANK=${OMPI_COMM_WORLD_RANK:-$SLURM_PROCID}

case $LRANK in
0|1) export UCX_NET_DEVICES=mlx5_0:1 ;;
2|3) export UCX_NET_DEVICES=mlx5_1:1 ;;
esac

echo "Rank $GRANK local $LRANK on $(hostname) UCX_NET_DEVICES=${UCX_NET_DEVICES:-auto}"

# "$@" keeps each argument intact
exec "$@"
