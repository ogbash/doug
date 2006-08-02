
# DOUG_UTILS_LIB_UNQUOTE(var_name,libs)
# "-labc -ldef" => "abc def"
AC_DEFUN([DOUG_UTILS_LIB_UNQUOTE],
[
_doug_utils_lib_qt=$2
_doug_utils_lib_nq=''
for l in $_doug_utils_lib_qt; do
    [_nq=`expr $l : '-l\([^ ]*\)'`]
    if test -z $_doug_utils_lib_nq; then
        _doug_utils_lib_nq=$_nq
    else
        _doug_utils_lib_nq="$_doug_utils_lib_nq $_nq"
    fi
done
$1=$_doug_utils_lib_nq
])

# DOUG_UTILS_LIB_QUOTE(var_name,libs)
# "abc def" => "-labc -ldef"
#AC_DEFUN([DOUG_UTILS_LIB_QUOTE],
#[
#
#])

# DOUG_UTILS_HEAD_TAIL(head_var,tail_var,string)
AC_DEFUN([DOUG_UTILS_HEAD_TAIL],
[
  _doug_utils_head_tail=$3
  set x $_doug_utils_head_tail; shift
  $1=$[1]; shift
  $2=$[*]
])
