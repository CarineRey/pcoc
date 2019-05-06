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
RUN pip install pandas

ENV LD_LIBRARY_PATH=/usr/local/lib/

# install bpp:

ENV bpp_core_version=aad2b8d
#405cab5 -> #aad2b8d

WORKDIR $bpp_dir/sources/bpp-core
RUN git clone https://github.com/BioPP/bpp-core . &&\
    git checkout $bpp_core_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

ENV bpp_seq_version=bb21df5
#32d9c67 -> #bb21df5

WORKDIR $bpp_dir/sources/bpp-seq
RUN git clone https://github.com/BioPP/bpp-seq . && \
    git checkout $bpp_seq_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install


ENV bpp_popgen_version=b788fd9
#77d712e -> #b788fd9

WORKDIR $bpp_dir/sources/bpp-popgen
RUN git clone https://github.com/BioPP/bpp-popgen . &&\
    git checkout $bpp_popgen_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install


ENV bpp_phyl_version=6b419d1 
#78d0235 -> #6b419d1 

WORKDIR $bpp_dir/sources/bpp-phyl
RUN git clone  https://github.com/BioPP/bpp-phyl . &&\
    git checkout $bpp_phyl_version &&\
    cmake  . && \
    (make -j 4 || make)  && \
    make install


ENV bppsuite_version=be18fc9
#77ccc0a -> #be18fc9

WORKDIR $bpp_dir/sources/bppsuite
RUN git clone https://github.com/BioPP/bppsuite . &&\
    git checkout $bppsuite_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

RUN pip install scipy
WORKDIR /data
COPY data/* /data/
COPY README.md /usr/local/etc/
COPY etc/entrypoint.sh /usr/local/bin/entrypoint.sh

CMD ["cat", "/usr/local/etc/README.md"]
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

COPY src/* /usr/local/bin/
