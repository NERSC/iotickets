#include <mpi.h>
#include <hdf5.h>
#include <cassert>
#include <math.h>
#include <vector>
#include <string>

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
  m_dsetDims[0] = m_numberOfArraysToWrite;
  m_dsetDims[1] = ARRAYDIM0;
  m_dsetDims[2] = ARRAYDIM1;
  m_dsetDims[3] = ARRAYDIM2;
  m_elemPerArray = ARRAYDIM0*ARRAYDIM1*ARRAYDIM2;
}

hid_t CreateAndFillFile::createFile() {
  hid_t file_id;
  hid_t	plist_id;
  MPI_Info info = MPI_INFO_NULL;  // could be way to give hints to hdf5
  assert((plist_id = H5Pcreate(H5P_FILE_ACCESS)) != -1);
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
  std::vector<int16_t> countArrays(count[0]*m_elemPerArray);

  hid_t memspace_id, filespace_id, plist_id;

  double t0 = MPI_Wtime();
  assert (0 <= (plist_id = H5Pcreate(H5P_DATASET_XFER)));
  assert (0 <= H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT));

  assert (0 <= (memspace_id = H5Screate_simple(DSETRANK, count, NULL)));
  assert (0 <= (filespace_id = H5Dget_space(dset_id)));

  assert (0 <= H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL));
  assert (0 <= H5Dwrite(dset_id, 
                        H5T_NATIVE_INT16, 
                        memspace_id, 
                        filespace_id, 
                        plist_id, &countArrays[0]));
  arraysWritten = count[0];
  bytesWritten = count[0]*sizeof(int16_t)*m_elemPerArray;
  assert (0 <= H5Pclose(plist_id));
  assert (0 <= H5Sclose(filespace_id));
  assert (0 <= H5Sclose(memspace_id));
  double totalTime = MPI_Wtime()-t0;
  double rateMB = bytesWritten/(1024*1024.0*totalTime);
  printf("rank=%d wrote %.2f MB block. time=%.2f sec rate=%.2f MB/sec\n", 
         rank, bytesWritten/(1024*1024.0), totalTime, rateMB);
}

bool CreateAndFillFile::run() {
  hid_t file_id = createFile();
  hid_t dset_id = createDataset(file_id);
  double t0 = MPI_Wtime();
  size_t bytesWritten;
  int arraysWritten;
  rankWrite(dset_id, m_mpi_rank, bytesWritten, arraysWritten);
  double thisRankFillTime = MPI_Wtime()-t0;
  printf("rank=%d worldsize=%d: run() fill of %d items %.5f seconds rate=%.2f MB/sec\n",
         m_mpi_rank, m_mpi_size, arraysWritten, thisRankFillTime,
         (bytesWritten/(1024*1024*thisRankFillTime)));

  assert (H5Dclose(dset_id)>=0);
  assert (H5Fclose(file_id)>=0);
  return true;
}

const char *getOutputFile() {
  // for hopper
  //return "/scratch/scratchdirs/davidsch/my_stripe_medium/out.h5";
  // for cori
  return "/global/cscratch1/sd/davidsch/my_stripe_medium/out.h5";
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int resultlen;
  MPI_Get_processor_name(hostname, &resultlen);
  const char *outputfile = getOutputFile();
  double gbToWrite = 4.0;
  double  t0 = MPI_Wtime();
  CreateAndFillFile createAndFillFile(outputfile, 
                                      gbToWrite, 
                                      comm, 
                                      hostname,
                                      mpi_rank, 
                                      mpi_size);

  bool result = createAndFillFile.run();
  double writeTime = MPI_Wtime()-t0;
  t0 = MPI_Wtime();
  createAndFillFile.check();
  double checkTime = MPI_Wtime()-t0;
  printf("%s: rank=%d worldsize=%d: wrote=%.2f GB sec=%.2f rate=%.2f MB/sec read=%.2f\n",
         hostname, mpi_rank, mpi_size, gbToWrite, writeTime, (gbToWrite*1024)/writeTime,
         checkTime);
  
    
  MPI_Finalize();
  return 0;
}

