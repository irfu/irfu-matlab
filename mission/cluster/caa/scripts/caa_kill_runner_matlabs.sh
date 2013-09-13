#!/bin/sh
#
# Kill active MATLAB sessions on the CAA Runners
# Useful for freeing up Matlab licenses or reloading ns_ops.

ssh caa 'sudo pkill -9 MATLAB'
ssh caa1 'sudo pkill -9 MATLAB'
ssh sipadan 'sudo pkill -9 MATLAB'

