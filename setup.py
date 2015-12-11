from numpy.distutils.core import Extension, setup

lz95_sources = ["fortran/lz95_model.f95", "fortran/lz95_model_tlm.f95",
    "fortran/lz95_model_adj.f95", "fortran/lz95_forecast.f95",
    "fortran/lz95_forecast_tlm.f95", "fortran/lz95_forecast_adj.f95"]

lz96_sources = ["fortran/lz96_model.f95", "fortran/lz96_model_tlm.f95",
    "fortran/lz96_model_adj.f95", "fortran/lz96_forecast.f95",
    "fortran/lz96_forecast_tlm.f95", "fortran/lz96_forecast_adj.f95"]

setup(name="lz95_fortran", ext_modules=[Extension(name="lz95_fortran",
    sources=lz95_sources)])

setup(name="lz96_fortran", ext_modules=[Extension(name="lz96_fortran",
    sources=lz96_sources)])
