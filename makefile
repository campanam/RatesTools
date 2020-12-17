ifndef INSTALL
	INSTALLDIR = $(realpath ${HOME}/ratestools)
else 
	INSTALLDIR = $(realpath $(INSTALL))
	ifeq ($(strip $(INSTALLDIR)),)
		INSTALLDIR = ${PWD}/$(INSTALL)
	endif
endif

BINDIR = bin

test:
	echo "RatesTools ruby executables will be installed into: $(INSTALLDIR)"
	if [[ ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo "$(INSTALLDIR) is in PATH."; else echo "$(INSTALLDIR) is not in PATH."; fi 

install: 
	if [ ! -d $(INSTALLDIR) ]; then mkdir $(INSTALLDIR); fi
	chmod +x $(BINDIR)/*.rb
	mv $(BINDIR)/*.rb $(INSTALLDIR)/
	mv ratestools.nf $(INSTALLDIR)/
	if [ ! -d ${HOME}/.profile ]; then mkfile -n 0 ${HOME}/.profile}; fi
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo 'export PATH="${PATH}:$(INSTALLDIR)"' >> ${HOME}/.profile; fi
	source "${HOME}/.profile"
	
.PHONY: test install

.SILENT: test install