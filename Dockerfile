FROM rust:latest

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --allow-unauthenticated --no-install-recommends \
        liblapack-dev \
        libopenblas-dev \
        gfortran \
        build-essential \
        sudo \
        npm

# Install rustfmt
RUN rustup component add rustfmt

# Make sure latest version of npm is installed
RUN npm install npm@latest -g

# Install rust wasm
RUN curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh

# Install cargo generate
RUN cargo install cargo-generate

# Copy the project inside the image
COPY . /root/mpc-rs

WORKDIR /root/mpc-rs
