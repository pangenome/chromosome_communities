#!/bin/bash

f=$1
o=$(mktemp louvain.XXXXX)

pafnet $f -e >$o.edgelist
pafnet $f -l >$o.names
louvain $o.edgelist $o.louvain >/dev/null

join <(tr ' ' '\t' <$o.louvain | sort -k 1b,1 ) \
     <(tr ' ' '\t' <$o.names | sort -k 1b,1) \
    | grep 'chm13\|grch38' | awk '{ print $2, $3 }' | grep -v _ \
    | sort -n | awk '{ c=$1 ; if (c == p) { x=x" "$2 } else { print x; x=$2 }; p = c; }'

rm -f $o.paf $o.edgelist $o.names $o.louvain $o
