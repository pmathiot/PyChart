#!/bin/bash
if [[ $# -ne 1 ]]; then echo 'E R R O R, usage: read_arg_list [input file]; exit'; exit; fi
sed -e :a -e N -e 's/\n/ /' -e ta $1
