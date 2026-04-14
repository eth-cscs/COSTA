FROM ubuntu:24.04 as builder

ENV DEBIAN_FRONTEND noninteractive

ENV FORCE_UNSAFE_CONFIGURE 1

ENV PATH="/spack/bin:${PATH}"

ENV MPICH_VERSION=3.4.3


RUN apt-get -y update

RUN apt-get install -y apt-utils

# install basic tools
RUN apt-get install -y --no-install-recommends gcc g++ gfortran git make unzip file \
  vim wget pkg-config python3-pip python3-dev cython3 python3-pythran curl tcl m4 cpio automake \
  xz-utils patch patchelf apt-transport-https ca-certificates gnupg software-properties-common perl tar bzip2 \
  liblzma-dev libbz2-dev

# get latest version of spack
RUN git clone -b releases/v1.1 https://github.com/spack/spack.git

COPY spack_repo /spack_repo

RUN spack repo add /spack_repo/costa

# set the location of packages built by spack
RUN spack config add config:install_tree:root:/opt/local

# find all external packages
RUN spack external find --all --exclude meson --exclude python

# find compilers
RUN spack compiler find

# install yq (utility to manipulate the yaml files)
RUN wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_arm64 && chmod a+x /usr/local/bin/yq

# install MPICH
RUN spack install mpich@${MPICH_VERSION} %gcc

# for the MPI hook
RUN echo $(spack find --format='{prefix.lib}' mpich) > /etc/ld.so.conf.d/mpich.conf
RUN ldconfig

# create environments for several configurations and install dependencies
RUN spack env create -d /COSTA-env && \
    spack -e /COSTA-env add "costa@master %gcc build_type=RelWithDebInfo +scalapack +tests +shared ^openblas threads=openmp" && \
    spack -e /COSTA-env develop -p /src costa@master && \
    spack -e /COSTA-env install --only=dependencies --fail-fast
