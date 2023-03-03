import pathlib
import platform
import setuptools
import setuptools.command.build_py
import setuptools.extension

dirname = pathlib.Path(__file__).resolve().parent

with open(dirname / "README.md") as file:
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
        for path in sorted((dirname / "astrometry.net" / "gsl-an" / "cblas").iterdir())
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


class BuildCommand(setuptools.command.build_py.build_py):
    def run(self):
        with open(
            dirname
            / "astrometry.net"
            / "include"
            / "astrometry"
            / "os-features-config.h",
            "w",
        ) as output:
            output.write(
                "\n".join(
                    (
                        "#define HAVE_NETPBM 0",
                        '#define AN_GIT_REVISION "0.89-7-g47849f04"',
                        '#define AN_GIT_DATE "Fri Jan 28 16:25:01 2022 -0500"',
                        '#define AN_GIT_URL "https://github.com/dstndstn/astrometry.net"\n',
                    )
                )
            )
        with open(dirname / "astrometry.net" / "gsl-an" / "config.h", "w") as output:
            output.write(
                "\n".join(
                    (
                        "#define GSL_DISABLE_DEPRECATED 1",
                        "#define HAVE_DECL_ACOSH 1",
                        "#define HAVE_DECL_ASINH 1",
                        "#define HAVE_DECL_ATANH 1",
                        "#define HAVE_DECL_EXPM1 1",
                        "#define HAVE_DECL_FEENABLEEXCEPT 0",
                        "#define HAVE_DECL_FESETTRAPENABLE 0",
                        "#define HAVE_DECL_FINITE 1",
                        "#define HAVE_DECL_FREXP 1",
                        "#define HAVE_DECL_HYPOT 1",
                        "#define HAVE_DECL_ISFINITE 1",
                        "#define HAVE_DECL_ISINF 1",
                        "#define HAVE_DECL_ISNAN 1",
                        "#define HAVE_DECL_LDEXP 1",
                        "#define HAVE_DECL_LOG1P 1",
                        "#define HAVE_DLFCN_H 1",
                        "#define HAVE_EXIT_SUCCESS_AND_FAILURE 1",
                        "#define HAVE_EXTENDED_PRECISION_REGISTERS 1",
                        "#define HAVE_INTTYPES_H 1",
                        "#define HAVE_LIBM 1",
                        "#define HAVE_MEMCPY 1",
                        "#define HAVE_MEMMOVE 1",
                        "#define HAVE_MEMORY_H 1",
                        "#define HAVE_PRINTF_LONGDOUBLE 1",
                        "#define HAVE_STDINT_H 1",
                        "#define HAVE_STDLIB_H 1",
                        "#define HAVE_STRDUP 1",
                        "#define HAVE_STRINGS_H 1",
                        "#define HAVE_STRING_H 1",
                        "#define HAVE_STRTOL 1",
                        "#define HAVE_STRTOUL 1",
                        "#define HAVE_SYS_STAT_H 1",
                        "#define HAVE_SYS_TYPES_H 1",
                        "#define HAVE_UNISTD_H 1",
                        "#define HAVE_VPRINTF 1",
                        '#define LT_OBJDIR ".libs/"',
                        '#define PACKAGE "gsl"',
                        '#define PACKAGE_BUGREPORT ""',
                        '#define PACKAGE_NAME "gsl"',
                        '#define PACKAGE_STRING "gsl 1.11"',
                        '#define PACKAGE_TARNAME "gsl"',
                        '#define PACKAGE_VERSION "1.11"',
                        "#define RELEASED",
                        "#define STDC_HEADERS 1",
                        '#define VERSION "1.11"',
                        "#if !HAVE_EXIT_SUCCESS_AND_FAILURE",
                        "#define EXIT_SUCCESS 0",
                        "#define EXIT_FAILURE 1",
                        "#endif",
                        "#if HAVE_EXTENDED_PRECISION_REGISTERS",
                        "#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))",
                        "#else",
                        "#define GSL_COERCE_DBL(x) (x)",
                        "#endif",
                        "#if !HAVE_DECL_HYPOT",
                        "#define hypot gsl_hypot",
                        "#endif",
                        "#if !HAVE_DECL_LOG1P",
                        "#define log1p gsl_log1p",
                        "#endif",
                        "#if !HAVE_DECL_EXPM1",
                        "#define expm1 gsl_expm1",
                        "#endif",
                        "#if !HAVE_DECL_ACOSH",
                        "#define acosh gsl_acosh",
                        "#endif",
                        "#if !HAVE_DECL_ASINH",
                        "#define asinh gsl_asinh",
                        "#endif",
                        "#if !HAVE_DECL_ATANH",
                        "#define atanh gsl_atanh",
                        "#endif",
                        "#if !HAVE_DECL_LDEXP",
                        "#define ldexp gsl_ldexp",
                        "#endif",
                        "#if !HAVE_DECL_FREXP",
                        "#define frexp gsl_frexp",
                        "#endif",
                        "#if !HAVE_DECL_ISINF",
                        "#define isinf gsl_isinf",
                        "#endif",
                        "#if !HAVE_DECL_ISFINITE",
                        "#define isfinite gsl_finite",
                        "#endif",
                        "#if !HAVE_DECL_FINITE",
                        "#define finite gsl_finite",
                        "#endif",
                        "#if !HAVE_DECL_ISNAN",
                        "#define isnan gsl_isnan",
                        "#endif",
                        "#ifdef __GNUC__",
                        "#define DISCARD_POINTER(p) do { ; } while(p ? 0 : 0);",
                        "#else",
                        "#define DISCARD_POINTER(p)",
                        "#endif",
                        "#if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)",
                        "#define GSL_RANGE_CHECK 0",
                        "#endif",
                        "#define RETURN_IF_NULL(x) if (!x) { return ; }\n",
                    )
                )
            )
        super().run()


setuptools.setup(
    name="astrometry",
    version="4.1.1",
    url="https://github.com/neuromorphicsystems/astrometry",
    author="ICNS, Alexandre Marcireau",
    author_email="alexandre.marcireau@gmail.com",
    description="Astrometry.net solver interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=["requests"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    cmdclass={
        "build_py": BuildCommand,
    },
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
