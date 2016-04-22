from setuptools import setup
from distutils.core import setup,Extension
from Cython.Build import cythonize

setup(
    name='SVelter',
    description='svelter',
    author='Xuefang Zhao',
    url='https://github.com/mills-lab/svelter',
    download_url='https://github.com/mills-lab/svelter',
    author_email='xuefzhao@umich.edu',
    version='1.1.0',
    packages=['svelter_sv'],
    scripts=["svelter_sv/svelter.py", "svelter_sv/SVelter.py"],
    ext_modules=cythonize("svelter_sv/*.pyx"),
    package_data={
        "SVelter": [
            "SVelter1.NullModel.Figure.a.r",
            "SVelter1.NullModel.Figure.b.r",
            "SVelter1.NullModel.Figure.c.r",
        ],
    },
)
#setup(**config)
