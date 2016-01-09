sudo apt-get -y build-dep r-base;
sudo add-apt-repository -y ppa:ubuntu-lxc/buildd-backports;
sudo apt-get update;
sudo apt-get install -y curl libcurl4-openssl-dev;
# zlib
wget http://zlib.net/zlib-1.2.8.tar.gz;
tar zxf zlib-1.2.8.tar.gz;
cd zlib-1.2.8;
./configure;
make;
sudo make install;
cd ../
  rm -rf zlib-1.2.8;
rm zlib-1.2.8.tar.gz;
wget --no-check-certificate https://stat.ethz.ch/R/daily/R-devel.tar.gz;
tar zxf R-devel.tar.gz;
cd R-devel;
./configure --with-x=no --without-system-pcre;
make;
sudo make install;
cd ../
rm -rf R-devel;
rm R-devel.tar.gz;
Rscript -e 'install.packages("devtools", dependencies=c("Depends", "Imports", "Suggests")';
Rscript -e 'install_github("biogram", dependencies=c("Depends", "Imports", "Suggests")';