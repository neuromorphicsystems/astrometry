import pathlib
import platform
import setuptools
import setuptools.extension

with open("README.md") as file:
    long_description = file.read()


extra_compile_args = []
extra_link_args = []
define_macros = []
include_dirs = [
    "astrometry.net/gsl-an",
    "astrometry.net/include",
    "astrometry.net/include/astrometry",
    "astrometry.net/util",
]
sources = [
    # gsl
    "astrometry.net/gsl-an/blas/blas.c",
    "astrometry.net/gsl-an/block/init.c",
    *(
        str(path.as_posix())
        for path in sorted(
            (pathlib.Path("astrometry.net") / "gsl-an" / "cblas").iterdir()
        )
        if path.suffix == ".c" and path.name != "hypot.c"
    ),
    "astrometry.net/gsl-an/err/error.c",
    "astrometry.net/gsl-an/err/stream.c",
    "astrometry.net/gsl-an/err/strerror.c",
    "astrometry.net/gsl-an/linalg/bidiag.c",
    "astrometry.net/gsl-an/linalg/householder.c",
    "astrometry.net/gsl-an/linalg/lu.c",
    "astrometry.net/gsl-an/linalg/qr.c",
    "astrometry.net/gsl-an/linalg/svd.c",
    "astrometry.net/gsl-an/matrix/copy.c",
    "astrometry.net/gsl-an/matrix/init.c",
    "astrometry.net/gsl-an/matrix/matrix.c",
    "astrometry.net/gsl-an/matrix/rowcol.c",
    "astrometry.net/gsl-an/matrix/submatrix.c",
    "astrometry.net/gsl-an/matrix/swap.c",
    "astrometry.net/gsl-an/matrix/view.c",
    "astrometry.net/gsl-an/permutation/init.c",
    "astrometry.net/gsl-an/permutation/permutation.c",
    "astrometry.net/gsl-an/permutation/permute.c",
    "astrometry.net/gsl-an/sys/coerce.c",
    "astrometry.net/gsl-an/sys/fdiv.c",
    "astrometry.net/gsl-an/sys/infnan.c",
    "astrometry.net/gsl-an/vector/copy.c",
    "astrometry.net/gsl-an/vector/init.c",
    "astrometry.net/gsl-an/vector/oper.c",
    "astrometry.net/gsl-an/vector/subvector.c",
    "astrometry.net/gsl-an/vector/swap.c",
    "astrometry.net/gsl-an/vector/vector.c",
    # libkd
    "astrometry.net/libkd/kdint_ddd.c",
    "astrometry.net/libkd/kdint_dds.c",
    "astrometry.net/libkd/kdint_ddu.c",
    "astrometry.net/libkd/kdint_dss.c",
    "astrometry.net/libkd/kdint_duu.c",
    "astrometry.net/libkd/kdint_fff.c",
    "astrometry.net/libkd/kdint_lll.c",
    "astrometry.net/libkd/kdtree.c",
    "astrometry.net/libkd/kdtree_dim.c",
    "astrometry.net/libkd/kdtree_fits_io.c",
    # qfits
    "astrometry.net/qfits-an/anqfits.c",
    "astrometry.net/qfits-an/qfits_byteswap.c",
    "astrometry.net/qfits-an/qfits_card.c",
    "astrometry.net/qfits-an/qfits_convert.c",
    "astrometry.net/qfits-an/qfits_error.c",
    "astrometry.net/qfits-an/qfits_float.c",
    "astrometry.net/qfits-an/qfits_header.c",
    "astrometry.net/qfits-an/qfits_memory.c",
    "astrometry.net/qfits-an/qfits_rw.c",
    "astrometry.net/qfits-an/qfits_table.c",
    "astrometry.net/qfits-an/qfits_time.c",
    "astrometry.net/qfits-an/qfits_tools.c",
    # solver
    "astrometry.net/solver/quad-utils.c",
    "astrometry.net/solver/solver.c",
    "astrometry.net/solver/tweak2.c",
    "astrometry.net/solver/verify.c",
    # util
    "astrometry.net/util/an-endian.c",
    "astrometry.net/util/bl.c",
    "astrometry.net/util/codekd.c",
    "astrometry.net/util/datalog.c",
    "astrometry.net/util/errors.c",
    "astrometry.net/util/fit-wcs.c",
    "astrometry.net/util/fitsbin.c",
    "astrometry.net/util/fitsfile.c",
    "astrometry.net/util/fitsioutils.c",
    "astrometry.net/util/fitstable.c",
    "astrometry.net/util/gslutils.c",
    "astrometry.net/util/healpix.c",
    "astrometry.net/util/index.c",
    "astrometry.net/util/ioutils.c",
    "astrometry.net/util/log.c",
    "astrometry.net/util/matchobj.c",
    "astrometry.net/util/mathutil.c",
    "astrometry.net/util/permutedsort.c",
    "astrometry.net/util/quadfile.c",
    "astrometry.net/util/sip.c",
    "astrometry.net/util/sip-utils.c",
    "astrometry.net/util/starkd.c",
    "astrometry.net/util/starutil.c",
    "astrometry.net/util/starxy.c",
    "astrometry.net/util/tic.c",
    # extension
    "astrometry_extension/astrometry_extension.c",
]

compiler = platform.python_compiler()
if compiler.startswith("Clang") or compiler.startswith("GCC"):
    extra_compile_args += ["-Wno-sign-compare"]

setuptools.setup(
    name="astrometry",
    version="2.0.0",
    url="https://github.com/neuromorphicsystems/astrometry",
    author="ICNS, Alexandre Marcireau",
    author_email="alexandre.marcireau@gmail.com",
    description="Astrometry.net solver interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    ext_modules=[
        setuptools.extension.Extension(
            "astrometry_extension",
            language="c",
            sources=sources,
            include_dirs=include_dirs,
            libraries=[],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            define_macros=define_macros,
        ),
    ],
    packages=["astrometry"],
)
