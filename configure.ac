AC_INIT([dmimgproject2], [4.13.5], [glotfeltyk@si.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AM_MAINTAINER_MODE([disable])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX


# -- Readline -- there has got to be a better way to do this ?

AC_CHECK_FUNC(tgetent,
   [READLINE_SUPPLIB=],
     [AC_CHECK_LIB(tinfo, tgetent,
     [READLINE_SUPPLIB=-ltinfo],
       [AC_CHECK_LIB(curses, tgetent,
         [READLINE_SUPPLIB=-lcurses],
         [AC_CHECK_LIB(ncurses, tgetent,
            [READLINE_SUPPLIB=-lncurses],
            [AC_CHECK_LIB(termcap, tgetent,
               [READLINE_SUPPLIB=-ltermcap],
               [AC_MSG_WARN(No termcap nor curses library found)])
            ])
         ])
      ])
   ])

AC_CHECK_LIB(readline,readline,[READLINELIB="-lreadline -lhistory -lncurses"],
	AC_MSG_ERROR([
You need the GNU readline library(ftp://ftp.gnu.org/gnu/readline/ ) to build
this program!
]), [$READLINE_SUPPLIB])

LIBS="$LIBS $PROCSTATLIB $READLINELIB $READLINE_SUPPLIB"

# ---------- CIAO PKG_CONFIG stuff
PKG_CONFIG_PATH=$prefix/lib/pkgconfig:$prefix/ots/lib/pkgconfig:$PKG_CONFIG_PATH
PKG_CONFIG_PATH=./pkgconfig:${PKG_CONFIG_PATH}
export PKG_CONFIG_PATH
PKG_CHECK_MODULES( [CIAO], [
  dmimgio
  dmfilters
  region
  ascdm
  cxcparam
  err
  ds
])


# ------------

AM_CONDITIONAL([LINUX], [test x`uname` == xLinux])


# -------------- Which files to configure ----------

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  wrapper/Makefile
  test/Makefile
])

AC_CONFIG_FILES([wrapper/dmimgproject2:wrapper/wrapper.sh],[chmod +x wrapper/dmimgproject2])

# -- Go! -----------

AC_OUTPUT
  
