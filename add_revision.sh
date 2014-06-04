#!/bin/sh

rev=`git rev-parse HEAD`
d=`git log -1 --format=%cd`
date="Last source update on: ${d} "


template=include/basic/global.template.h
global=include/basic/global.h
eval "sed -e 's/SVN_VERSION_TAG/$rev/g' -e 's/SVN_COMMIT_DATE_TAG/$date/g' $template " > $global
echo "${rev}"
echo "${date}"
