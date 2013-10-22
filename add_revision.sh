#!/bin/sh

rev=`svn info | grep Revision`
date=`svn info | grep Date | awk '{day=(0+$8);print "Last source update on " $9" " day " " $10 " at " $5 }' | sed -e s/\)//g`


template=include/basic/global.template.h
global=include/basic/global.h
eval "sed -e 's/SVN_VERSION_TAG/$rev/g' -e 's/SVN_COMMIT_DATE_TAG/$date/g' $template " > $global
echo "${rev}"
