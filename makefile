#----------------------------------------------------------------------------------------
# Michael G. Campana and Ellie E. Armstrong, 2020-2022
# Smithsonian Institution and Stanford University

# RatesTools Makefile 0.5.3

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2022. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. *bioRxiv*.
# doi: 10.1101/2022.07.18.500472.
#----------------------------------------------------------------------------------------

ifndef INSTALL
	INSTALLDIR = $(realpath ${HOME})/ratestools
else 
	INSTALLDIR = $(realpath $(INSTALL))
	ifeq ($(strip $(INSTALLDIR)),)
		INSTALLDIR = ${PWD}/$(INSTALL)
	endif
endif

BINDIR = bin

test:
	echo "RatesTools ruby executables will be installed into: $(INSTALLDIR)"
	if [[ ":${PATH}:" == *":$(INSTALLDIR)"* ]]; then echo "$(INSTALLDIR) is in PATH."; else echo "$(INSTALLDIR) is not in PATH."; fi 

install: 
	if [ ! -d $(INSTALLDIR) ]; then mkdir $(INSTALLDIR); fi
	chmod +x $(BINDIR)/*.rb
	chmod +x $(BINDIR)/*.R
	chmod +x $(BINDIR)/*.sh
	mv $(BINDIR)/*.rb $(INSTALLDIR)/
	mv $(BINDIR)/*.R $(INSTALLDIR)/
	mv $(BINDIR)/*.sh $(INSTALLDIR)/
	chmod +x ratestools.nf
	mv ratestools.nf $(INSTALLDIR)/
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)"* ]]; then \
		if [[ -f ${HOME}/.bash_profile ]]; then \
			if [[ `grep 'export PATH="$$PATH' ${HOME}/.bash_profile | wc -l` -gt 0 ]]; then \
				sed -i '' 's#export PATH="$$PATH#export PATH="$$PATH:$(INSTALLDIR)#' ${HOME}/.bash_profile; \
			else echo 'export PATH="$$PATH:$(INSTALLDIR)"' >> ${HOME}/.bash_profile; \
			fi; \
		else \
			echo 'export PATH="$$PATH:$(INSTALLDIR)"' >> ${HOME}/.bash_profile; \
		fi; \
		if [[ -f ${HOME}/.bashrc ]]; then \
			if [[ `grep 'export PATH="$$PATH' ${HOME}/.bashrc | wc -l` -gt 0 ]]; then \
				sed -i '' 's#export PATH="$$PATH#export PATH="$$PATH:$(INSTALLDIR)#' ${HOME}/.bashrc; \
			else echo 'export PATH="$$PATH:$(INSTALLDIR)"' >> ${HOME}/.bashrc; \
			fi; \
		else \
			echo 'export PATH="$$PATH:$(INSTALLDIR)"' >> ${HOME}/.bashrc; \
		fi; \
		if [[ -f ${HOME}/.zshrc ]]; then \
			if [[ `grep 'export PATH="$$PATH' ${HOME}/.zshrc | wc -l` -gt 0 ]]; then \
				sed -i '' 's#export PATH="$$PATH#export PATH="$$PATH:$(INSTALLDIR)#'  ${HOME}/.zshrc; \
			else echo 'export PATH="$$PATH:$(INSTALLDIR)"' >> ${HOME}/.zshrc; \
			fi; \
		else \
			echo 'export PATH="$$PATH:$(INSTALLDIR)"' >> ${HOME}/.zshrc; \
		fi; \
	fi
	
.PHONY: test install

.SILENT: test install
