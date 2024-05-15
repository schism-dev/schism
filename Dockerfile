# SPDX-FileCopyrightText: 2022-2023 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de

FROM phusion/baseimage:jammy-1.0.1

LABEL description="SCHISM Docker environment based on Ubuntu"
LABEL author="Carsten Lemmen <carsten.lemmen@hereon.de>"
LABEL license="CC0-1.0"
LABEL copyright="2022-2023 Helmholtz-Zentrum hereon GmbH"

# Arguments can be passed via the --build-arg key=value command to the 
# docker build command.  The default values are set below to mpich, oldio=on,
# and superbee 
ARG COMMUNICATOR="mpich"
ARG OLDIO="ON"
ARG TVD_LIM="SB"

RUN apt update && apt -qy install cmake gcc-11 python3 \
    python-is-python3 lib${COMMUNICATOR}-dev libmetis-dev libnetcdf-dev \
    libnetcdff-dev libparmetis-dev git
ENV PATH="/usr/lib64/${COMMUNICATOR}/bin:${PATH}"

WORKDIR /usr/src

RUN git clone  --branch master --depth 1 https://github.com/schism-dev/schism /usr/src/schism
RUN mkdir -p /usr/src/schism/build
RUN cmake -S /usr/src/schism/src -B /usr/src/schism/build -DOLDIO=${OLDIO} -DTVD_LIM=${TVD_LIM}
RUN make -C /usr/src/schism/build -j8 pschism
