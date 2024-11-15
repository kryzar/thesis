# Support repository for my PhD thesis

This is the support repository for my [PhD thesis](thesis.pdf). We explain:
- [how to get the SageMath version we use in our thesis](#running-sagemath);
- [how to run our benchmarks](#running-benchmarks).

## Running SageMath

We explain how to run the SageMath version of the thesis.

### Run SageMath from the already built Docker Hub image

The docker image associated to [the Dockerfile](Dockerfile) has already been built, and was pushed to [the Docker Hub](https://hub.docker.com/repository/docker/kryzar/thesis-sagemath). It can directly be used. First of all, *pull* the image:

```bash
docker pull kryzar/thesis-sagemath
```

To run it:

```bash
docker run --rm --interactive kryzar/thesis-sagemath
```

> [!WARNING]
> The image is about 10 GB in size.

### Run your own Docker image

To build the Docker image associated to the Dockerfile, clone
the repository:

```bash
git clone https://github.com/kryzar/thesis
```

then go to the newly created directory (which contains the Dockerfile):

```bash
cd thesis
```

and build the image using `docker build` (the name `leudiere-thesis-sagemath` can be freely changed):

```bash
docker build -t leudiere-thesis-sagemath .
```

> [!NOTE]
> It is possible to build the image with `docker-buildx`, [an extension of Docker](https://docs.docker.com/reference/cli/docker/buildx/), running `docker buildx build -t leudiere-thesis-sagemath .`. This may offer better performances. The extension must be installed; on Debian and Ubuntu, use `sudo apt install docker-buildx`.

Once the image is built, run it with `docker run`; the following opens a SageMath shell:

```bash
docker run --rm --interactive kryzar/thesis-sagemath
```

Optionally, to push the image to the Docker Hub, run:

```bash
docker tag leudiere-thesis-sagemath:latest kryzar/thesis-sagemath:latest
docker login -u kryzar  # Will prompt for a password; use an access token
docker push kryzar/thesis-sagemath:latest
```

### Manually build SageMath on your computer

One can avoid using Docker and manually build SageMath from source. The compilation instructions can be found [here](https://doc.sagemath.org/html/en/installation/source.html) and [here](https://github.com/sagemath/sage/?tab=readme-ov-file#instructions-to-build-from-source). We also mention that the [SageMath Developer Guide](https://doc.sagemath.org/html/en/developer/index.html) is valuable learning resource.

> [!IMPORTANT]
> Be sure to have all the dependencies installed. Their list for each platform (Linux, macOS, Windows) is available [here](https://github.com/sagemath/sage/?tab=readme-ov-file#instructions-to-build-from-source). However, we stress that those prerequisites are heavy and numerous. In order to keep the userspace clean, using Docker could be recommended.

> [!IMPORTANT]
> Be sure to use [the branch `thesis` of `kryzar/sagemath`](https://github.com/kryzar/sage/tree/thesis) and not [the branch `develop` of `sagemath/sage`](https://github.com/sagemath/sage) (default branch to build SageMath).

Roughly, the sequence of commands is the following:

> [!WARNING]
> Run those at your own risk.

```bash
git clone -c core.symlinks=true --filter blob:none \
          --origin origin --branch thesis --tags   \
          https://github.com/kryzar/sage.git sage
cd sage
make configure
./configure
make build -j64  # Parallelize on 64 processors
```
Once this has finished, run:

```bash
./sage
```

> [!TIP]
> As per [the installation guide](https://github.com/sagemath/sage), we highly recommend parallelizing the build process. This can be achieved by running `make build -jN` (where `N` is the number of processors); for example, `make build -j64` builds on 64 parallel processors, while `make build -j8` builds on 8 parallels processors.

> [!NOTE]
> The command `make build` does not compile the whole SageMath suite, but only the main executable. For instance, the documentation is not built. This saves significant time, as well as space.

## Running benchmarks

We give the code for two classes of benchmarks:
- The benchmarks of Chapters 4 and 5. The code is available in [benchmarks/chapters-4-5](benchmarks/chapters-4-5). To reproduce our tests, run the module `benchmarks.sage` with an appropriate version of SageMath. The results as well as the console logs are written to various files; the global variable `WORKDIR` (line 17) should therefore be set to an appropriate value. We have fixed the random seed so that the exact same random isogenies and endomorphisms are used.
- The benchmarks of Chapter 6. The code is available in [benchmarks/chapter-6](benchmarks/chapter-6). These were originally written by [Xavier Caruso](https://xavier.caruso.ovh/) ([@xcaruso](https://github.com/xcaruso)) for the Software Presentation [*Drinfeld Modules in SageMath*](https://dl.acm.org/doi/10.1145/3614408.3614417) (ACM Communications in Computer Algebra, Volume 57, Issue 2, 2023).

