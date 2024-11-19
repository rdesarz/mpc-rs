FROM rust:latest

ARG C_USER
ARG C_UID
ARG C_GID

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --allow-unauthenticated --no-install-recommends \
        liblapack-dev
        gfortran
        build-essential


RUN groupadd -g $C_GID $C_USER && \
    useradd -m -d /home/$C_USER -g $C_GID -s /bin/bash -u $C_UID $C_USER && \
    adduser $C_USER sudo

USER $C_USER
WORKDIR /home/$C_USER 
