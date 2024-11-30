FROM rust:latest

ARG C_USER
ARG C_UID
ARG C_GID

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --allow-unauthenticated --no-install-recommends \
        liblapack-dev \
        gfortran \
        build-essential \
        npm

RUN groupadd -g $C_GID $C_USER && \
    useradd -m -d /home/$C_USER -g $C_GID -s /bin/bash -u $C_UID $C_USER && \
    adduser $C_USER sudo

COPY --chown=$C_UID:$C_GID . /home/$C_USER/mpc-rs

USER $C_USER
WORKDIR /home/$C_USER/mpc-rs
