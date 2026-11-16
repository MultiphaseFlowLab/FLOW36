#!/bin/bash
# Vega GPU nodes: InfiniBand card per GPU not set yet (UCX picks the devices itself).
# TODO: set UCX_NET_DEVICES per local rank once the node topology is known
# (nvidia-smi topo -m; ls /sys/class/infiniband), as in Leonardo/binder_leo.sh
LRANK=$(( ${OMPI_COMM_WORLD_LOCAL_RANK} ))

echo "Rank $OMPI_COMM_WORLD_RANK local $LRANK on $(hostname) UCX_NET_DEVICES=${UCX_NET_DEVICES:-auto}"

# "$@" keeps each argument intact
exec "$@"
