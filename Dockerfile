# syntax=docker/dockerfile:1

# Builds a statically linked (musl) siblinggwas binary on Alpine and ships it in an empty image.
#
#   docker build -t mrcieu/siblinggwas .
#   docker run --rm -v "$PWD:/data" -w /data mrcieu/siblinggwas siblinggwas --raw ... --out ...

FROM rust:1.94-alpine AS build
RUN apk add --no-cache musl-dev

WORKDIR /build
COPY Cargo.toml Cargo.lock NOTICE ./
COPY src src
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    cargo build --release --locked \
    && strip target/release/siblinggwas

FROM scratch
COPY --from=build /build/target/release/siblinggwas /usr/local/bin/siblinggwas
COPY LICENSE NOTICE /usr/share/doc/siblinggwas/
ENV PATH=/usr/local/bin
CMD ["siblinggwas", "--help"]
