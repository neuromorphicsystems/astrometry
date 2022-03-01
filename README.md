# astrometry
A Python wrapper for core Astrometry.net functions

```sh
rm -rf astrometry.egg-info; rm -rf astrometry_extension.cpython-39-darwin.so; rm -rf build; python3 setup.py develop
```

Pre-build
```sh
python3 prebuild.py
```

```sh
clang-format -i astrometry_extension/astrometry_extension.c; CC="ccache clang" python3 setup.py develop --disable-lto;
```

- output_logodds_threshold: immediately accept matches whose log-odds are larger than this value
- tune_up_logodds_threshold: tune-up matches whose log-odds are larger than this value and accept them if their post-tune-up log-odds are larger than output_logodds_threshold

In both cases, accepted matches are tuned up again. Because log-odds are compared with the thresholds before this last tune-up, the final log-odds are often significantly larger than `output_logodds_threshold`.

Set `tune_up_logodds_threshold` to a value larger than or equal to `output_logodds_threshold` to disable the first tune-up.

See below for a description of the tuning up algorithm (from *astrometry.net/include/astrometry/tweak2.h*) and `tweak2` in *astrometry.net/solver/tweak2.c* for the implementation.

Given an initial WCS solution, compute SIP polynomial distortions using an annealing-like strategy. That is, it finds matches between image and reference catalog by searching within a radius, and that radius is small near a region of confidence, and grows as you move away.  That makes it possible to pick up more distant matches, but they are downweighted in the fit. The annealing process reduces the slope of the growth of the matching radius with respect to the distance from the region of confidence.


series_4100: 355.54 MB
series_4200: 33.78 GB
series_5000: 76.24 GB
series_5200: 36.14 GB
series_6000: 1.20 GB
series_6100: 1.58 GB


## Contribute

To format the code, run:
```sh
clang-format -i astrometry_extension/astrometry_extension.c
```

For local development, run:
```sh
python3 setup.py develop
# use CC="ccache clang" python3 setup.py develop
# to speed up incremental builds
```

## Publish

1. Bump the version number in *setup.py*.

2. Install Cubuzoa in a different directory (https://github.com/neuromorphicsystems/cubuzoa) to build pre-compiled versions for all major operating systems. Cubuzoa depends on VirtualBox (with its extension pack) and requires about 75 GB of free disk space.
```
cd cubuzoa
python3 cubuzoa.py provision
python3 cubuzoa.py build /path/to/astrometry
```

3. Install twine
```
pip3 install twine
```

4. Upload the compiled wheels and the source code to PyPI:
```
python3 setup.py sdist --dist-dir wheels
python3 -m twine upload wheels/*
```
