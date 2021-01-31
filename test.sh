#!/bin/bash

make schwinger2peD
cd jobs

cd test
rm -rf {ckpoint,topo_charge,pion_corr,vacuum_trace,plots}
./schwinger2peD.sh
