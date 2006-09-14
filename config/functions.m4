dnl -*- shell-script -*-
dnl
dnl This file is taken from LAM/MPI distribution
dnl
dnl Copyright (c) 2001-2004 The Trustees of Indiana University.  
dnl                         All rights reserved.
dnl Copyright (c) 1998-2001 University of Notre Dame. 
dnl                         All rights reserved.
dnl Copyright (c) 1994-1998 The Ohio State University.  
dnl                         All rights reserved.
dnl 
dnl This file is part of the LAM/MPI software package.  For license
dnl information, see the LICENSE file in the top level directory of the
dnl LAM/MPI source distribution.
dnl
dnl $Id: lam_functions.m4,v 1.17 2004/01/20 03:41:47 jsquyres Exp $
dnl

AC_DEFUN([LOG_MSG],[
# 1 is the message
# 2 is whether to put a prefix or not
if test -n "$2"; then
    echo "configure:__oline__: $1" >&5
else
    echo $1 >&5
fi])dnl

dnl #######################################################################

AC_DEFUN([LOG_FILE],[
# 1 is the filename
if test -n "$1" -a -f "$1"; then
    cat $1 >&5
fi])dnl

dnl #######################################################################

AC_DEFUN([LOG_COMMAND],[
# 1 is the command
# 2 is actions to do if success
# 3 is actions to do if fail
echo "configure:__oline__: $1" >&5
$1 1>&5 2>&1
lam_status=$?
LOG_MSG([\$? = $lam_status], 1)
if test "$lam_status" = "0"; then
    unset lam_status
    $2
else
    unset lam_status
    $3
fi])dnl
