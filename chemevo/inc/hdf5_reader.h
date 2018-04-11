#ifndef HDF5
#define HDF5
#include <H5Cpp.h>
using namespace H5;
//=============================================================================
// HDF5 writers
//=============================================================================
template<typename T> struct get_hdf5_data_type
{   static H5::PredType type()
    {
        //static_assert(false, "Unknown HDF5 data type");
        return H5::PredType::NATIVE_DOUBLE;
    }
};
template<> struct get_hdf5_data_type<double>{
    FloatType type{H5::PredType::NATIVE_DOUBLE};};
/**
 * @brief Read 1D vector from HDF5
 *
 * @param fin input HDF5 file
 * @param data 1D vector
 * @param name name of data
 */

template<class c>
void hdf5_read_1D_vector(H5File &fin, std::vector<c> &data, std::string name){
    DataSet dataset = fin.openDataSet(name);
    DataSpace dataspace = dataset.getSpace();
    H5T_class_t dataClass = dataset.getTypeClass();
    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    data.resize(dims_out[0]);
    dataset.read(data.data(),PredType::NATIVE_DOUBLE,dataspace);
}
/**
 * @brief Write 1D vector to HDF5
 *
 * @param fout output HDF5 file
 * @param data 1D vector
 * @param name name of data
 */

template<class c>
void hdf5_write_1D_vector(H5File &fout, std::vector<c> &data, std::string name){
    hsize_t dim(1);
    DataSpace dspace(1,&dim);
    get_hdf5_data_type<c> typ;
    VarLenType dtype(&typ.type);
    DataSet dataset = fout.createDataSet(name,dtype,dspace);
    hvl_t vl[dim];
    vl[0].len = data.size();
    vl[0].p = &data[0];
    dataset.write(vl,dtype);
}
/**
 * @brief Write 2D vector of vectors to HDF5
 *
 * @param fout output HDF5 file
 * @param data 2D vector of vectors
 * @param name name of data
 */
template<class c>
void hdf5_write_2D(H5File &fout,std::vector<std::vector<c>> &data, std::string name){
    hsize_t dim(data.size());
    DataSpace dspace(1,&dim);
    get_hdf5_data_type<c> typ;
    VarLenType dtype(&typ.type);
    DataSet dataset = fout.createDataSet(name,dtype,dspace);
    hvl_t vl[dim];
    for(hsize_t i=0;i<dim;++i){
        vl[i].len = data[i].size();
        vl[i].p = &data[i][0];
    }
    dataset.write(vl,dtype);
    return;
}
#endif
