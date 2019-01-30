#!/bin/bash

_=`which sw_vers`
if [ 0 -eq $? ] # macOS
then
	set -e
	name=`sw_vers -productName`
else
	set -e
	name=`uname -s`
fi
echo "${name}" | awk '{ gsub(" ", "") ; print tolower($0) }'
