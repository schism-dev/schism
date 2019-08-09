#!/bin/sh

#
# This script is used to make a new release of GOTM.
# The release can be both a 'stable' and 'devel' release.
# The script should not be executed directly - but via the 
# make file targets <devel|stable>

# BUGS:
# the syncronisation with the gotm.net web page should be done
# via rsync instead of using scp and ssh
#

#
# $Id
#

[ "$USER" = "kbk" ] || { echo "Only kbk can make new releases" ; exit 1; }

release_type=$1
release_version=$2

release_name=gotm-$release_version
base_dir=/public/ftp/pub/gotm-releases
release_dir=$base_dir/$release_type
tarfile=$release_name.tar.gz

TAG=v`echo $release_version | tr . _`
BRANCH=$TAG

RHOST=bolding@gotm.net
RDIR=/data/kamel/domains/gotm.net/src/

export CVSROOT=$USER@gate:/public/cvs
export CVS_RSH=ssh

if [ -d $release_dir/$release_name ] ; then

   echo
   echo $release_name" has already been released"
   echo "update VERSION in Makefile"
   echo
   exit 1

fi

if [ "$release_type" = "stable" ] ; then
   cvs tag $TAG
   CVS2CL="cvs2cl -b -F v3_2_0 --no-ancestors"
fi

if [ "$release_type" = "devel" ] ; then
   cvs tag $TAG
   CVS2CL="cvs2cl -F trunk"
fi

if [ "$release_type" = "branch" ] ; then
   cvs tag -b $TAG
   CVS2CL="cvs2cl -b -F $BRANCH --no-ancestors"
fi


$CVS2CL && mkdir -p $release_dir/$release_name/include/ && mv ChangeLog VERSION $release_dir/$release_name && mv include/version.h $release_dir/$release_name/include/

cd $release_dir && cvs export -r $TAG -d $release_name gotm && tar -cvzf $tarfile $release_name/ && ln -sf $release_name.tar.gz gotm_$release_type.tar.gz && ln -sf $release_name/Changelog

#rsync -av --include=*.tar.gz --exclude=* $base_dir/ $RHOST:$RDIR

scp -p $release_dir/$tarfile  $release_dir/Changelog $RHOST:$RDIR/$release_type/
ssh $RHOST \( cd $RDIR/$release_type \; rm gotm-$release_type.tar.gz \; ln -s $tarfile gotm-$release_type.tar.gz \)

exit 0
