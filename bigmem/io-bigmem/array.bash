#!/usr/bin/env bash
set -e;

id=$[$SLURM_ARRAY_TASK_ID - 0];
printf -v id "_%04d" $id;
cd /global/cscratch1/sd/jialin/io-ticket/nerscio/io-bigmem/${id};
bash ${id}.bash;

