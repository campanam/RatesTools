# RatesTools  

**Michael G. Campana & Ellie E. Armstrong, 2019-2020**  
Smithsonian Conservation Biology Institute  
Stanford University  

Pipeline to calculate de novo mutation rates from parent-offspring trios  

## License  
This software is available under  

## Installation and Configuration  
### Install Nextflow and Ruby  
RatesTools requires [Nextflow](https://www.nextflow.io/) and [Ruby](http://www.ruby-lang.org). Basic instructions for installing these languages are copied below. We recommend installing Ruby using the [Ruby Version Manager](https://rvm.io). See the official language documentation should you need help installing these languages.  

Install Nextflow: `wget -qO- https://get.nextflow.io | bash`  
Install the latest Ruby using Ruby Version Manager: `curl -sSL https://get.rvm.io | bash -s stable --ruby`  

### Install the RatesTools Package  
Clone the repository: `git clone https://github.com/campanam/RatesTools`  
Install the scripts: `cd RatesTools; make install`  
*By default, RatesTools scripts will be installed into the ~/ratestools directory. If you wish to change the default directory specify the INSTALL parameter, e.g.:* `make INSTALL=/path/to/some/dir install`  

