#!/bin/bash

IFS='\n'
TAGS=`git tag -l --points-at HEAD`
HASH=`git rev-parse --short --verify HEAD`
VERSION="DEV"
if [ -n "${TAGS}" ]
then
	for name in ${TAGS}
	do
		if [[ "${name}" == v* ]]
		then
			VERSION="${name:1}"
			break
		fi
	done
fi
echo ${VERSION} ${HASH}
