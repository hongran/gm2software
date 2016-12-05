#!/bin/bash

#cd $(dirname $(readlink -f $0))/odb.d
cd ./odb.d

for script in `ls [0-9]?-*.cmd`
do
    odbedit -e TestMagnet -c @${script}
done
