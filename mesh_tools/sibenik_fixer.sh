#!/bin/bash

#This fixes everything except the quads
grep '^[fv] ' sibenik.obj | tr -s ' ' | sed 's/\/[0-9]\+//g;s/ $//g' | python3 triangler.py > fixed_sibenik.obj
