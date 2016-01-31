SUBDIRS = utls curvlinear curvidet

all:
	@for f in $(SUBDIRS) ; do ( cd $$f; $(MAKE) all) done

install:
	@for f in $(SUBDIRS) ; do ( cd $$f; $(MAKE) install) done

clean:
	@for f in $(SUBDIRS) ; do ( cd $$f; $(MAKE) clean) done
	rm -rf include
	rm -rf lib/LNX86
