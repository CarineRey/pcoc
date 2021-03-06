FROM debian:stretch
MAINTAINER Carine Rey carine.rey@ens-lyon.org

RUN apt-get update && \
    apt-get install --no-install-recommends -qy \
       git \
       make \
       cmake \
       python2.7-minimal \
       ca-certificates \
       curl \
       g++\
       gcc \
       wget \
       gnupg2 \
       dirmngr \
       python-pip \
       python-numpy \
       python-pandas \
       python-setuptools \
       python-dev \
       sudo \
       less \
       xvfb \
       gosu \
       pyqt4-dev-tools \
       xauth \
       libcurl4-openssl-dev \
       libxml2-dev \
       libssl-dev \
       libcairo2-dev

RUN pip install ete3==3.0.0b35
RUN pip install biopython

ENV LD_LIBRARY_PATH=/usr/local/lib/

# # Install R
# RUN echo """deb http://cloud.r-project.org/bin/linux/debian stretch-cran34/""" >> /etc/apt/sources.list \
#     && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys FCAE2A0E115C3D8A
# RUN apt-get update && apt-get install --no-install-recommends -qy r-base r-base-dev \
#     && rm -rf /var/lib/apt/lists/* \
#     && echo 'install.packages(c("ggplot2", "plyr", "reshape2", "readr", "cowplot"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
#     && Rscript /tmp/packages.R

# install bpp:
WORKDIR $bpp_dir/sources/bpp-core
RUN git clone https://github.com/BioPP/bpp-core . &&\
    git checkout 405cab5 &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

WORKDIR $bpp_dir/sources/bpp-seq
RUN git clone https://github.com/BioPP/bpp-seq . && \
    git checkout 32d9c67 &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

WORKDIR $bpp_dir/sources/bpp-popgen
RUN git clone https://github.com/BioPP/bpp-popgen . &&\
    git checkout 77d712e &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

WORKDIR $bpp_dir/sources/bpp-phyl
RUN git clone  https://github.com/BioPP/bpp-phyl . &&\
    git checkout  78d0235 &&\
    cmake  . && \
    (make -j 4 || make)  && \
    make install

WORKDIR $bpp_dir/sources/bppsuite
RUN git clone https://github.com/BioPP/bppsuite . &&\
    git checkout 77ccc0a &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install


COPY data/* /data/
COPY README.md /usr/local/etc/
COPY etc/entrypoint.sh /usr/local/bin/entrypoint.sh

CMD ["cat", "/usr/local/etc/README.md"]
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

COPY src/* /usr/local/bin/
