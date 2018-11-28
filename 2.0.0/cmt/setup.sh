# echo "setup SYJ_Lambda_c 2.0.0 in /workfs/bes/suyj/BES3_WorkArea/7.0.3"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtSYJ_Lambda_ctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtSYJ_Lambda_ctempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=SYJ_Lambda_c -version=2.0.0 -path=/workfs/bes/suyj/BES3_WorkArea/7.0.3  -no_cleanup $* >${cmtSYJ_Lambda_ctempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=SYJ_Lambda_c -version=2.0.0 -path=/workfs/bes/suyj/BES3_WorkArea/7.0.3  -no_cleanup $* >${cmtSYJ_Lambda_ctempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtSYJ_Lambda_ctempfile}
  unset cmtSYJ_Lambda_ctempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtSYJ_Lambda_ctempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtSYJ_Lambda_ctempfile}
unset cmtSYJ_Lambda_ctempfile
return $cmtsetupstatus

