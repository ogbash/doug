
# empty for now, add SCRIPTS= later
dist_noinst_SCRIPTS=amdeps.py manager.py scripts.py \
	gridplot.py autotools.py compareSolutions.py dougtestmysql.py dougtest.py \
	dougtesttar.py doxscript.py doxygen.py svnscripts.py testscript.py

EXTRA_DIST=doug modules

dist-hook:
	rm -rf `find $(distdir)/doug -name .svn`
	rm -rf `find $(distdir)/modules -name .svn`
