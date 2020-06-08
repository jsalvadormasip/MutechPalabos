#ifndef HDF_WRAPPER
#define HDF_WRAPPER

#ifdef HDF5

#include "hdf5.h"
#include "core/globalDefs.h"
#include "stdlib.h"
#include "string.h"
#include <vector>


#define DATASETNAME 	"IntArray" 
#define NX     8                      /* dataset dimensions */
#define NY     6 
#define RANK   1


#pragma GCC diagnostic ignored "-Wunused-function"
//variable global
static hsize_t file_id_g = 0;

static void write_string_hdf5(hid_t file_id, std::string const& string, int mpi_rank) {

	hid_t       dataset_id, dataspace_id;  /* identifiers */
	hsize_t		nb = string.size();
	herr_t      status;
	/* Open an existing file. */
	
	dataspace_id = H5Screate_simple(1, &nb, NULL);
	dataset_id   = H5Dcreate2(file_id, "/meta_data", H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if (mpi_rank == 0) {
		status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, string.c_str());
		if (status < 0) {
			printf("Writing string with hdf5 didn't work, error code: %d \n", status);
		}
	}
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	
}

//todo make it parallel.. maybe
static char *read_string_hdf5() {

	hid_t       dset_id;  /* identifiers */
	herr_t      status;
	dset_id = H5Dopen2(file_id_g, "/meta_data", H5P_DEFAULT);
	hid_t dspace = H5Dget_space(dset_id);
	hsize_t dim, maxdims;
	H5Sget_simple_extent_dims(dspace, &dim, &maxdims);

	char *xm_read = (char*)malloc((int)maxdims + 1);
	status = H5Dread(dset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, xm_read);
	if (status < 0) {
		printf("Reading string with hdf5 didn't work, error code: %d \n", status);
	}
	xm_read[maxdims] = 0;
	H5Dclose(dset_id);
	return xm_read;
}

template<typename T = char >
static void write_parallel_hdf5(const char *file_name, const char *data_set_path, std::vector<plb::plint> const& my_block_id, std::vector<plb::plint> const& i_offset, std::vector<std::vector<char> >&data, std::string const& metadata, int mpi_rank, MPI_Comm comm, bool opened = false){
	
	hid_t h5_type = H5T_C_S1, h5_element_size = sizeof(T);

	if (typeid(T) == typeid(float)) {
		h5_type = H5T_NATIVE_FLOAT;
		//h5_element_size = sizeof(float);
	}
	if (typeid(T) == typeid(double)) {
		h5_type = H5T_NATIVE_DOUBLE;
		printf("h5_type %d \n", h5_type);
		//h5_element_size = sizeof(double);
	}
	
    hid_t   file_id, dset_id, plist_id;         /* file and dataset identifiers */
    hid_t   filespace, memspace;      /* file and memory dataspace identifiers */
	herr_t	status;
	hsize_t	data_piece_size;
	hsize_t	offset; 
	hsize_t	data_set_size = i_offset.back()/h5_element_size;
	
	if (my_block_id.size() == 0) {
		data_piece_size = 0;
		offset = 0;

	} else if ( my_block_id[0] > 0){
		data_piece_size  = i_offset[my_block_id.back()] - i_offset[my_block_id[0]-1];
		offset  	     = i_offset[my_block_id[0]-1]/h5_element_size;
	}else{
		data_piece_size  = i_offset[my_block_id.back()];
		offset  	     = 0;
	}
	//int		offset_output	 = i_offset + piece_size;

	//printf(" offset %d | data_piece_size %d data_set_size %d rand %d |my_block_id[0] %d my_block_id[n_blocks-1] %d \n", (int)offset, (int)data_piece_size, (int)data_set_size, mpi_rank, my_block_id[0], my_block_id[n_blocks-1]);
	
	char *data_total = (char*)malloc((int)data_piece_size);
	data_piece_size /= h5_element_size;
	
	int idx_acc = 0;
	
	for (int i = 0; i < (int)my_block_id.size(); i++) {
		
		int block_size;
		
		if(my_block_id[i] > 0){
			block_size = i_offset[my_block_id[i]] - i_offset[my_block_id[i]-1];
		}else{
			block_size = i_offset[0];
		}
		
		//printf("block_size %d idx_acc %d | rank %d data_piece_size %d \n", block_size, idx_acc, mpi_rank, data_piece_size);
		
		for(int j = 0; j < block_size; j++){
			data_total[idx_acc + j] = data[i][j];
		}
		idx_acc += block_size;
	}
	
	MPI_Info info  = MPI_INFO_NULL;
	//creat the file
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
	//std::cout << "file name " << file_name << std::endl;
	if(opened){
		file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
		
	}else{
		file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	}

    H5Pclose(plist_id);

	//std::cout << "metadata size | " << metadata << " | ran " << mpi_rank << std::endl;
	if(metadata[0])write_string_hdf5(file_id, metadata, mpi_rank);
	
	//creat dataset   
	filespace = H5Screate_simple(1, &data_set_size, NULL);
    memspace  = H5Screate_simple(1, &data_piece_size, NULL);

    dset_id = H5Dcreate(file_id, data_set_path, h5_type, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
	
	//creat pieces to be written in parallel
	filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &data_piece_size, NULL);
	
	//write data in parallel
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite(dset_id, h5_type, memspace, filespace, plist_id, data_total);


	if (status < 0) {
		printf("hdf5 write didnt work, error code: %d \n", status);
		return;
	}
	free(data_total);
    //Close/release resources.
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
	
    H5Pclose(plist_id);
	H5Fclose(file_id);
}

static void open_hdf_file(const char *file_name, MPI_Comm comm) {
	hid_t   plist_id;
	MPI_Info info = MPI_INFO_NULL;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, comm, info);
	file_id_g = H5Fopen(file_name, H5F_ACC_RDONLY, plist_id);
	H5Pclose(plist_id);
}

static void close_hdf_file() {
	H5Fclose(file_id_g);
}

//for this to work the file has to be already open (the file descriptor is a global var), 
//this is an optimization it allows to read data and metadata without having to close and reopen the file
static std::vector< std::vector<char> > read_parallel_hdf5(std::vector<plb::plint> const& my_block_id, std::vector<plb::plint> const& i_offset, int mpi_rank, MPI_Comm comm) {
	
	hid_t   dset_id, plist_id;         // file and dataset identifiers //
	hid_t   filespace, memspace;      // file and memory dataspace identifiers //
	herr_t	status;
	hsize_t	data_piece_size;
	hsize_t	offset;
	//hsize_t	data_set_size = i_offset.back();

	//std::cout << " my_block_Id " << my_block_id.size() << " " << mpi_rank << std::endl;

	if( my_block_id.size() == 0){
		data_piece_size = 0;
		offset = 0;

	} else if  (my_block_id[0] > 0) {
		data_piece_size = i_offset[my_block_id.back()] - i_offset[my_block_id[0] - 1];
		offset = i_offset[my_block_id[0] - 1];

	} else {
		data_piece_size = i_offset[my_block_id.back()];
		offset = 0;
	}

	//MPI_Info info = MPI_INFO_NULL;

	/////////////////////allocat memory for the data to be read////
	char *data_total = (char*)malloc((int)data_piece_size);
	
	//creat dataset   
    dset_id = H5Dopen2(file_id_g, "binary_blob", H5P_DEFAULT);
	
	filespace = H5Dget_space(dset_id);
	memspace  = H5Screate_simple(1, &data_piece_size, NULL);
	
	//creat pices to be read in parallel//
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &data_piece_size, NULL);
	
	//read data in parallel
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dread(dset_id, H5T_C_S1, memspace, filespace, plist_id, data_total);

	if (status < 0) {
		printf("hdf5 read didnt work, error code: %d \n", status);
		std::vector< std::vector<char> > empty;
		return empty;
	}
	//printf("data_piece_size %d \n", data_piece_size);
	
    //Close/release resources.
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

	std::vector< std::vector<char> >data(my_block_id.size());

	int idx_acc = 0;
	for (int i = 0; i < (int)my_block_id.size(); i++) {

		int block_size;

		if (my_block_id[i] > 0) {
			block_size = i_offset[my_block_id[i]] - i_offset[my_block_id[i] - 1];
		}
		else {
			block_size = i_offset[0];
		}
		data[i] = std::vector<char>(block_size);
		for (int j = 0; j < block_size; j++) {
			data[i][j] = data_total[idx_acc + j];
		}
		idx_acc += block_size;
	}
	free(data_total);
	return data;
}
#endif
#endif