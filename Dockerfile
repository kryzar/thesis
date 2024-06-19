FROM sagemath/sagemath-dev:develop

# Build instructions:
#     https://github.com/sagemath/sage/?tab=readme-ov-file#instructions-to-build-from-source

RUN sudo git clone -c core.symlinks=true --filter blob:none \
                   --origin origin --branch thesis --tags   \
                   https://github.com/kryzar/sage.git /sage

WORKDIR /sage

RUN sudo make configure
RUN sudo env
RUN sudo ./configure --enable-build-as-root
RUN sudo make -j64

ENTRYPOINT ["/sage/sage"]
