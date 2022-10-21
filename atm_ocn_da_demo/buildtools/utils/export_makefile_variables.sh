#!/bin/bash

EXPORT_MAKEFILE_VARIABLES_MYPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null && pwd )"

tmpfile=$(mktemp)
gmake -pn -f $EXPORT_MAKEFILE_VARIABLES_MYPATH/../Makefile.display_variables 2> /dev/null | grep -A1 "^# makefile"| grep -v "^#\|^--" | sort | uniq > $tmpfile

while read line
do
    var=$(echo "$line" | sed "s#^\([^[:space:]]*\) *\:*= *\(.*\)\$#\1#g")
    value=$(echo "$line" | sed "s#^\([^[:space:]]*\) *\:*= *\(.*\)\$#\2#g")
    if [ "$var" == ".DEFAULT_GOAL" ]; then continue; fi
    if [ "$var" == "MAKEFLAGS" ]; then continue; fi
    echo "export $var=\"$value\""
    eval "export $var=\"$value\""
done < "$tmpfile"

unlink $tmpfile
