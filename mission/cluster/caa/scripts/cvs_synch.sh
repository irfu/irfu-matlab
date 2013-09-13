#!/bin/sh
#
# Synchronize all of the CVS-controlled files required for CAA processing

  cd /home/chris/devel/irfu-matlab/caa/caa-control
  cvs up
  make all
  make txt
  cvs -f commit -m '' ns_ops.xml
  cvs -f commit -m '' ns_ops_c?.dat
  cvs -f commit -m '' manual_problems_c?.dat
  ssh hq 'cd /home/www/cluster/efw/ops; cvs up;'
  ssh db 'cd /data/cluster/irfu-matlab; cvs up; cd /data/cluster/caa-control/; cvs up'

