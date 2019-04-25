FROM debian:stretch
LABEL org.opencontainers.image.authors=carine.rey@ens-lyon.org

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
       libcairo2-dev \
       libeigen3-dev

RUN pip install ete3==3.0.0b35
RUN pip install biopython==1.76
RUN pip install pandas
RUN pip install scipy

ENV LD_LIBRARY_PATH=/usr/local/lib/

# install bpp:
RUN sudo apt-get install texinfo -y

ENV bpp_core_version=fa5da67
#405cab5 -> #aad2b8d -> #fa5da67

WORKDIR $bpp_dir/sources/bpp-core
RUN git clone https://github.com/BioPP/bpp-core . &&\
    git checkout $bpp_core_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

ENV bpp_seq_version=3930130
#32d9c67 -> #bb21df5 -> #3930130

WORKDIR $bpp_dir/sources/bpp-seq
RUN git clone https://github.com/BioPP/bpp-seq . && \
    git checkout $bpp_seq_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

ENV bpp_popgen_version=4bc260d
#77d712e -> #b788fd9 -> #4bc260d

WORKDIR $bpp_dir/sources/bpp-popgen
RUN git clone https://github.com/BioPP/bpp-popgen . &&\
    git checkout $bpp_popgen_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

ENV bpp_phyl_version=266d544
#78d0235 -> #4bbd6b1 -> #266d544

WORKDIR $bpp_dir/sources/bpp-phyl
RUN git clone  https://github.com/BioPP/bpp-phyl . &&\
    git checkout $bpp_phyl_version &&\
    cmake  . && \
    (make -j 4 || make)  && \
    make install

ENV bppsuite_version=a82cfc3
#77ccc0a -> #be18fc9 -> #a82cfc3

WORKDIR $bpp_dir/sources/bppsuite
RUN git clone https://github.com/BioPP/bppsuite . &&\
    git checkout $bppsuite_version &&\
    cmake  . && \
    (make -j 4 || make) && \
    make install

WORKDIR /data
COPY data/* /data/
COPY README.md /usr/local/etc/
COPY etc/entrypoint.sh /usr/local/bin/entrypoint.sh

COPY src/* /usr/local/bin/
CMD ["cat", "/usr/local/etc/README.md"]
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

COPY src/* /usr/local/bin/
