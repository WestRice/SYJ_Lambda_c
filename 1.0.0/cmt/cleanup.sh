# echo "cleanup SYJ_Lambda_c 1.0.0 in /workfs/bes/suyj/BES3_WorkArea/7.0.3"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtSYJ_Lambda_ctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtSYJ_Lambda_ctempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=SYJ_Lambda_c -version=1.0.0 -path=/workfs/bes/suyj/BES3_WorkArea/7.0.3  $* >${cmtSYJ_Lambda_ctempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=SYJ_Lambda_c -version=1.0.0 -path=/workfs/bes/suyj/BES3_WorkArea/7.0.3  $* >${cmtSYJ_Lambda_ctempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtSYJ_Lambda_ctempfile}
  unset cmtSYJ_Lambda_ctempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtSYJ_Lambda_ctempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtSYJ_Lambda_ctempfile}
unset cmtSYJ_Lambda_ctempfile
return $cmtcleanupstatus

