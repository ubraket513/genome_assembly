from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
except ImportError as exc:
    raise SystemExit(
        "Cython is required to build the optional accelerator. "
        "Install it with: python -m pip install '.[perf]'"
    ) from exc


extensions = [
    Extension(
        "genome_assembly._cython_backend",
        sources=["genome_assembly/_cython_backend.pyx"],
    )
]


setup(
    name="genome-assembly",
    packages=["genome_assembly"],
    ext_modules=cythonize(
        extensions,
        annotate=False,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        },
        nthreads=4,
    )
)
