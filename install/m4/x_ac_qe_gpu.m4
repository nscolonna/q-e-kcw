# Copyright (C) 2026 Quantum ESPRESSO Foundation
#

AC_DEFUN([X_AC_QE_GPU], [

AC_ARG_WITH(gpu,
   [AS_HELP_STRING([--with-gpu],
       [(cuda|omp5) Use "cuda" for NVidia, "omp5" for AMD (default:cuda)])],
   [if test "$withval" = "omp5" ; then
      with_gpu=2
   else
      with_gpu=1
   fi],
   [with_gpu=0])
 
AC_ARG_WITH(cuda,
   [AS_HELP_STRING([--with-cuda],
       [obsolete, use --with-gpu=cuda instead])],
   [if test "$withval" != "no" ; then
       AC_MSG_WARN([--with-cuda is obsolete, use --with-gpu=cuda instead])
       with_gpu=1
   fi],
   [with_gpu=0])
   
])
