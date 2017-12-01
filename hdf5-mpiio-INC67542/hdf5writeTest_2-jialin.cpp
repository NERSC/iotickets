#include <mpi.h>
#include "getopt.h"
#include <hdf5.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#define NAME_MAX 255
char cb_buffer_size[NAME_MAX];
char cb_nodes[NAME_MAX];
char filename[NAME_MAX];

class CreateAndFillFile {
public:
  CreateAndFillFile(const char *outputFileName,
                    double gbToWrite,
                    MPI_Comm &comm,
                    const char *hostname,
                    int mpi_rank,
                    int mpi_size);

  bool run();
  bool check() { return true; };

protected:
  enum {DSETRANK=4, ARRAYRANK=3, ARRAYDIM0=32, ARRAYDIM1=185, ARRAYDIM2=388};

  hid_t createFile();
  hid_t createDataset(hid_t file_id);
  void rankWrite(hid_t dset_id, int rank, size_t &bytesWritten, int &arraysWritten);
  hid_t openFile();
  hid_t openDataset(hid_t file_id);

private:
  const char *m_outputFileName;
  double m_gbToWrite;
  MPI_Comm &m_comm;
  std::string m_hostname;
  int m_mpi_rank;
  int m_mpi_size;
  int m_numberOfArraysToWrite;
  hsize_t m_dsetDims[DSETRANK];
  int m_elemPerArray;
  typedef int16_t Array[ARRAYDIM0][ARRAYDIM1][ARRAYDIM2];
};


CreateAndFillFile::CreateAndFillFile(const char *outputFileName,
                                     double gbToWrite,
                                     MPI_Comm &comm,
                                     const char *hostname,
                                     int mpi_rank,
                                     int mpi_size) :
  m_outputFileName(outputFileName)
  , m_gbToWrite(gbToWrite)
  , m_comm(comm)
  , m_hostname(hostname)
  , m_mpi_rank(mpi_rank)
  , m_mpi_size(mpi_size)
 {
  double gb2bytes = 1<<30;
  m_numberOfArraysToWrite = int(ceil((gbToWrite*gb2bytes)/double(sizeof(Array))));
  if(mpi_rank==0) printf("m_dsetDims: %d\n",m_numberOfArraysToWrite);
  m_dsetDims[0] = m_numberOfArraysToWrite;
  m_dsetDims[1] = ARRAYDIM0;
  m_dsetDims[2] = ARRAYDIM1;
  m_dsetDims[3] = ARRAYDIM2;
  m_elemPerArray = ARRAYDIM0*ARRAYDIM1*ARRAYDIM2;
}

hid_t CreateAndFillFile::createFile() {
  hid_t file_id;
  hid_t	plist_id;
  MPI_Info info;  // could be way to give hints to hdf5
  MPI_Info_create(&info);
  MPI_Info_set(info, "cb_buffer_size", cb_buffer_size);
  MPI_Info_set(info, "cb_nodes", cb_nodes);
  assert((plist_id = H5Pcreate(H5P_FILE_ACCESS)) != -1);
  //assert (H5Pset_fapl_mpiposix(plist_id, m_comm, info) >= 0); //only available in 1.8.12 and before
  assert (H5Pset_fapl_mpio(plist_id, m_comm, info) >= 0);
  assert((file_id = H5Fcreate(m_outputFileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id)) >= 0);
  assert (H5Pclose(plist_id) >=0);
  return file_id;
}

hid_t CreateAndFillFile::createDataset(hid_t file_id) {
  const char *dsetName = "data";
  hid_t filespace_id;
  hid_t dset_id;
  assert( 0 <= (filespace_id = H5Screate_simple(DSETRANK, m_dsetDims, NULL)) );
  assert( 0 <= (dset_id = H5Dcreate(file_id, dsetName, H5T_NATIVE_INT16, filespace_id,
                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)));
  assert (0 <= H5Sclose(filespace_id));
  return dset_id;
}

void CreateAndFillFile::rankWrite(hid_t dset_id, int rank, size_t &bytesWritten, int &arraysWritten) { 

  hsize_t offset[DSETRANK] = {rank * (m_numberOfArraysToWrite/m_mpi_size), 0, 0, 0};
  hsize_t count[DSETRANK] = {0,ARRAYDIM0,ARRAYDIM1,ARRAYDIM2};
  count[0] = std::min(hsize_t(m_numberOfArraysToWrite/m_mpi_size),
                      hsize_t(m_numberOfArraysToWrite-offset[0]));
  //printf("rank %d;offset:%d;count:%d\n",rank,offset[0],count[0]);
  std::vector<int16_t> countArrays(count[0]*m_elemPerArray);

  hid_t memspace_id, filespace_id, plist_id;

  //double t0 = MPI_Wtime();
  assert (0 <= (plist_id = H5Pcreate(H5P_DATASET_XFER)));
  assert (0 <= H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE));
  assert (0 <= (memspace_id = H5Screate_simple(DSETRANK, count, NULL)));
  assert (0 <= (filespace_id = H5Dget_space(dset_id)));

  assert (0 <= H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL));
  //double t0 = MPI_Wtime();
  assert (0 <= H5Dwrite(dset_id, 
                        H5T_NATIVE_INT16, 
                        memspace_id, 
                        filespace_id, 
                        plist_id, &countArrays[0]));
  arraysWritten = count[0];
  bytesWritten = count[0]*sizeof(int16_t)*m_elemPerArray;
  assert (0 <= H5Sclose(filespace_id));
  assert (0 <= H5Sclose(memspace_id));
}

bool CreateAndFillFile::run() {
  hid_t file_id = createFile();
  hid_t dset_id = createDataset(file_id);
  double t0 = MPI_Wtime();
  size_t bytesWritten;
  int arraysWritten;
  printf("rank=%d, starting at %.5f \n",m_mpi_rank,t0);
  rankWrite(dset_id, m_mpi_rank, bytesWritten, arraysWritten);
  double t1=MPI_Wtime();
  printf("rank=%d, stoping at %.5f \n",m_mpi_rank,t1);
  double thisRankFillTime = t1-t0;
  /*printf("rank=%d worldsize=%d: run() fill of %d items %.5f seconds rate=%.2f MB/sec\n",
         m_mpi_rank, m_mpi_size, arraysWritten, thisRankFillTime,
         (bytesWritten/(1024*1024*thisRankFillTime)));
  */
  assert (H5Dclose(dset_id)>=0);
  assert (H5Fclose(file_id)>=0);
  return true;
}

const char *getOutputFile() {
 return filename;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  double dsizeinGB=6.0;
  strncpy(cb_buffer_size,"16777216", NAME_MAX);
  strncpy(cb_nodes, "2", NAME_MAX);
  strncpy(filename, "./testoutput/out.h5", NAME_MAX);
  int c;
  while ((c = getopt (argc, argv, "f:b:n:x:")) != -1)
    switch (c)
      {
      case 'f':
        strncpy(filename, optarg, NAME_MAX);
        break;
      case 'b':
        strncpy(cb_buffer_size,optarg, NAME_MAX);
        break;
      case 'n':
        strncpy(cb_nodes, optarg, NAME_MAX);
        break;
      case 'x':
        dsizeinGB = strtoull(optarg, NULL, 10);
        break;
      default:
        break;
      }
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int resultlen;
  MPI_Get_processor_name(hostname, &resultlen);
  const char *outputfile = getOutputFile();
  double gbToWrite = dsizeinGB;
  double  t0 = MPI_Wtime();
  CreateAndFillFile createAndFillFile(outputfile, 
                                      gbToWrite, 
                                      comm, 
                                      hostname,
                                      mpi_rank, 
                                      mpi_size);

  bool result = createAndFillFile.run();
  double writeTime = MPI_Wtime()-t0;
  if(mpi_rank==0)
  printf("%s: rank=%d,nproc:%d, totalwrote=%.2f GB, sec=%.2f, rate=%.2f MB/sec",
         hostname, mpi_rank, mpi_size, gbToWrite, writeTime, (gbToWrite*1024)/writeTime);   
  MPI_Finalize();
  return 0;
}
