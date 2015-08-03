
MODLIST = lndpm lndcal lndsr lndsrbm 6sV-1.0B

all:
	@for i in $(MODLIST); do \
	(cd $$i; make; cd ..) \
	done

install:
	@for i in $(MODLIST); do \
	(cd $$i; make install; cd ..) \
	done

clean:
	@for i in $(MODLIST); do \
	(cd $$i; make clean; cd ..) \
	done
