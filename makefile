#----------------------------------------------------------------------------------------
# Michael G. Campana and Ellie E. Armstrong, 2020-2022
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Campana, M.G. & E.E. Armstrong. 2020. RatesTools: Pipeline to calculate de novo
# mutation rates from parent-offspring trios. Smithsonian Institution and Stanford
# University. <https://github.com/campanam/RatesTools>.
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
	if [[ ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo "$(INSTALLDIR) is in PATH."; else echo "$(INSTALLDIR) is not in PATH."; fi 

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
	if [ ! -d ${HOME}/.profile ]; then echo '' > ${HOME}/.profile}; fi
	if [ ! -d ${HOME}/.bash_profile ]; then echo '' > ${HOME}/.bash_profile}; fi
	if [ ! -d ${HOME}/.bashrc ]; then echo '' > ${HOME}/.bashrc}; fi
	if [ ! -d ${HOME}/.zprofile ]; then echo '' > ${HOME}/.zprofile}; fi
	if [ ! -d ${HOME}/.zshrc ]; then echo '' > ${HOME}/.zshrc}; fi
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo 'export PATH="${PATH}:$(INSTALLDIR)"' >> ${HOME}/.profile; fi
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo 'export PATH="${PATH}:$(INSTALLDIR)"' >> ${HOME}/.bash_profile; fi
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo 'export PATH="${PATH}:$(INSTALLDIR)"' >> ${HOME}/.bashrc; fi
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo 'export PATH="${PATH}:$(INSTALLDIR)"' >> ${HOME}/.zprofile; fi
	if [[ ! ":${PATH}:" == *":$(INSTALLDIR)":* ]]; then echo 'export PATH="${PATH}:$(INSTALLDIR)"' >> ${HOME}/.zshrc; fi
	source "${HOME}/.profile"
	source "${HOME}/.zprofile"
	source "${HOME}/.bash_profile"
	source "${HOME}/.bashrc"
	source "${HOME}/.zshrc"
	
.PHONY: test install

.SILENT: test install