# echo "cleanup SYJ_Lambda_c 2.0.0 in /workfs/bes/suyj/BES3_WorkArea/7.0.3"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtSYJ_Lambda_ctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtSYJ_Lambda_ctempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=SYJ_Lambda_c -version=2.0.0 -path=/workfs/bes/suyj/BES3_WorkArea/7.0.3  $* >${cmtSYJ_Lambda_ctempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=SYJ_Lambda_c -version=2.0.0 -path=/workfs/bes/suyj/BES3_WorkArea/7.0.3  $* >${cmtSYJ_Lambda_ctempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtSYJ_Lambda_ctempfile}
  unset cmtSYJ_Lambda_ctempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtSYJ_Lambda_ctempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtSYJ_Lambda_ctempfile}
unset cmtSYJ_Lambda_ctempfile
exit $cmtcleanupstatus

