#!/bin/bash

set -e

_=`which sw_vers`
if [ 0 -eq $? ] # macOS
then
	name=`sw_vers -productName`
else
	name=`uname -s`
fi
echo "${name}" | awk '{ gsub(" ", "") ; print tolower($0) }'
