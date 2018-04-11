#include "boost_python_interface.h"
#include "edf.h"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(moments_overload,sb15_edf::moments , 4, 5);

BOOST_PYTHON_MODULE(edf_py)
{
    docstring_options doc_options(true);

    // ========================================================================
    // EDF
    // ========================================================================
  class_<sb15_edf, boost::noncopyable>("edf",
            "EDF from Sanders & Binney (2015)\n"
            "\n"
            "All potential classes must have a Phi(x) and Forces(x) routine in C++. These are called by __call__ and Forces respectively in python",init<Potential_JS*,Actions_AxisymmetricFudge_InterpTables*,optional<VecDoub,VecDoub,bool,bool,int>>())
      .def("__call__", &sb15_edf::full_DF_Z,
           "Compute EDF ({x,y,z,vx,vy,vz},age,Z)\n"
            "\n"
            "Args:\n"
            "    param1: np.array of Cartesian position x.\n"
            "    param2: age tau.\n"
            "    param3: metallicity Z.\n"
            "\n"
            "Returns:\n"
            "    edf at point\n"
        "")
        .def("readParams",&sb15_edf::readParams,
             "Read parameters in file")
        .def("moments",&sb15_edf::moments,
             moments_overload(args("R","z","tau","Z"),
             "Compute density, means and dispersions at (R,z,age,metallicity)"));
}

