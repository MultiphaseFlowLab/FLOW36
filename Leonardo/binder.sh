#!/bin/bash
case $(( ${OMPI_COMM_WORLD_LOCAL_RANK} )) in
0) export UCX_NET_DEVICES=mlx5_0:1 ;;
1) export UCX_NET_DEVICES=mlx5_1:1 ;;
2) export UCX_NET_DEVICES=mlx5_2:1 ;;
3) export UCX_NET_DEVICES=mlx5_3:1 ;;
esac

echo Launching on $UCX_NET_DEVICES

$*

