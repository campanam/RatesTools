ifndef INSTALL
INSTALLDIR = $(realpath ${HOME}/ratestools)
else
INSTALLDIR =  $(realpath $(INSTALL))
endif

BINDIR = bin

test:
	echo "RatesTools ruby executables will be installed into: $(INSTALLDIR)"
	if [[ ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo "$(INSTALLDIR) is in PATH."; else echo "$(INSTALLDIR) is not in PATH."; fi 

install: 
	if [ ! -d $(INSTALLDIR) ]; then mkdir $(INSTALLDIR); fi
	chmod +x $(BINDIR)/*.rb
	mv $(BINDIR)/*.rb $(INSTALLDIR)/

.PHONY: test install

.SILENT: test install