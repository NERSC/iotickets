#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <fstream>

#define EIGEN_DEFAULT_TO_ROW_MAJOR

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "Lasso.h"
#include "Numbers.h"
#include "hdf5.h"

#define OUTPUTFILE	"output.h5"
#define DATASET_X 	"/data/X"
#define DATASET_Y	"/data/y"
#define MASTER		0 
#define edison		24
#define cori		32


using namespace Eigen;
using namespace std;

float is_not_NaN (float x) {if (!isnan(x)) return x; else return 0;}
pointer_to_unary_function <float,float> floorObject (floor) ;

int main (int argc, char *argv[]) 
{ 
   
   Lasso lasso; 
   Numbers num;

   
   char INPUTFILE[30]; 
   strcpy(INPUTFILE, argv[1]); 
   
   if (argc < 4) {cout << "Usage: ./bolbo_admm [input hdf5]  [--nMP] [--bootE] [--bootS]" << endl;}
   
   int world_rank, world_size;		//for world_comm rank and size of the process'
   int maxBoot;
   int nMP = atoi(argv[2]); 
   int nbootE = atoi(argv[3]); 
   int nbootS = atoi(argv[4]);
   float rndfrctL=0.8, rndfrct=0.8, cvlfrct=0.9;  
   int m,n; 			//Data size
   int m_frac, L_frac;
   int nrnd = 10; 
   int bgdOpt = 1;
   int seed = 1234;

   double start_sizeTime, end_sizeTime, start_commTime, end_commTime, start_loadTime, end_loadTime, start_winInitTime, start_distTime, end_distTime, start_compTime, end_compTime, start_las1Time, end_las1Time, start_bca1Time, end_bca1Time, start_las2Time, end_las2Time, start_bca2Time, end_bca2Time, start_olsTime, end_olsTime, start_bca3Time, end_bca3Time, start_saveTime, end_saveTime, total_startTime, total_endTime, start_bagTime, end_bagTime, start_preolsTime, end_preolsTime, start_LTTime, end_LTTime, start_ifTime, end_ifTime, start_elseTime, end_elseTime, start_sprtTime, end_sprtTime, start_permTime, end_permTime, end_setTime, start_setTime, start_perm2Time, end_perm2Time, start_preloopTime, end_preloopTime, end_prelas1Time, start_denseTime, end_denseTime, end_bogusTime, start_ccloopTime, end_ccloopTime; //time calculation

    /*MPI init and global communicator*/ 

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Group world_group; 
    MPI_Comm_group(MPI_COMM_WORLD, &world_group); 
    MPI_Status *mpi_status; 


    if (world_rank == MASTER)
   	total_startTime = MPI_Wtime();

    if (nbootE != nbootS) { maxBoot = nbootE > nbootS ? nbootE : nbootS; }
    else 		  { maxBoot = nbootE; }
 
    int processes = nMP * maxBoot; 	
 

    if (world_size < processes)
 	cout << "\nThe no. of processes do not match !!! Processes size should be atleast" << processes; 



    int nodes = ceil(world_size/edison);

    MPI_Comm comm_ADMM; 
    int color;

    //if (nodes == 1) { color = world_size; } else { color  = floor (world_rank/nodes);  }

   color = world_rank / edison; 

    cout << "color = " << color << endl;  

    MPI_Comm_split (MPI_COMM_WORLD, color, world_rank, &comm_ADMM); 
    int admm_rank, admm_size;
    MPI_Comm_rank(comm_ADMM, &admm_rank);
    MPI_Comm_size(comm_ADMM, &admm_size);

    MPI_Barrier (comm_ADMM); 

    cout << "admm RANK/SIZE= " << admm_rank << "/" << admm_size << endl; 


    VectorXi admm_root(nodes);
    int root = 0;	
   
    for (int i=0; i<nodes; i++) {
    	admm_root(i) = root;
	root+=edison; 
    }  

    MPI_Group admm_root_group; 
    MPI_Group_incl(world_group, nodes, admm_root.data(), &admm_root_group);
    MPI_Comm comm_ROOT; 
    //MPI_Comm_create_group(comm_ADMM, admm_root_group, 0, &comm_ROOT); 
    MPI_Comm_create_group(MPI_COMM_WORLD, admm_root_group, 0, &comm_ROOT);
    MPI_Status status;
    /*if (nodes != admm_size)
	cout << "\nRequested pool size and ADMM comm size do not match!!" << endl; */

    int root_rank = -1, root_size = -1; 

    if (MPI_COMM_NULL != comm_ROOT) {
  	 MPI_Comm_rank(comm_ROOT, &root_rank);
   	 MPI_Comm_size(comm_ROOT, &root_size); 
	 cout << "WORLD RANK/SIZE : " << world_rank << "/" << world_size << "---- ROOT RANK/SIZE : " << root_rank << "/" << root_size << endl;
    }
    

    VectorXi mp, mpid; 
    mp.setLinSpaced(nMP,0,nMP).colwise().replicate(nodes);
    sort(mp.data(), mp.data()+mp.size());
    mpid = mp.colwise().replicate(maxBoot);

    /*ofstream myfile ("data/mp_id.dat");
        if (myfile.is_open())
        {
         myfile << mpid;
         myfile.close();
        }*/
    
    int mp_id = mpid(admm_rank); 
    int bt_id = floor(world_rank/nMP);   
 
    if (world_rank == MASTER)
	start_sizeTime = MPI_Wtime(); 
  

    if (world_rank == MASTER) {
   	hid_t file_id, dataset_id;  
    	file_id = H5Fopen(INPUTFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    	dataset_id = H5Dopen(file_id, DATASET_X, H5P_DEFAULT);
	hid_t dspace = H5Dget_space(dataset_id);
        const int ndims = H5Sget_simple_extent_ndims(dspace);
	hsize_t dims[ndims];
    	H5Sget_simple_extent_dims(dspace, dims, NULL);
	m = dims[0];
   	n = dims[1]; 
        H5Dclose(dataset_id);
        H5Fclose(file_id); 
     }
   
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int m_node = floor(m/nodes); 

    int chunksize = floor (m/nodes);
    int offset = chunksize;

    if (world_rank == MASTER) {
 	end_sizeTime = MPI_Wtime() - start_sizeTime; 
	cout << "Data size: (" << m << "," << n << ")" << endl;
	cout << "\nData size broadcast: " << end_sizeTime << endl; 
    }

    VectorXi row_ids(m), rows_node(m_node), wheredata(m), wheredata_node(m_node), nodeid_gen(nodes);
    int node_id;   
  
    if (world_rank == MASTER) { 
	row_ids.setLinSpaced(m,0,m); 
	random_shuffle (row_ids.data(), row_ids.data()+row_ids.size());

	wheredata = row_ids / m_node;
		
	for (int i=0; i < row_ids.size(); i++)
		row_ids(i) = row_ids(i) % m_node; 

	if (nodes > 2)
		nodeid_gen.setLinSpaced(nodes,0,nodes);
	else if (nodes == 2)
		nodeid_gen << 0, 1; 
	else
		nodeid_gen << 0; 
   
	 /*ofstream myfile ("row_ids.dat");
    	if (myfile.is_open())
    	{
       	 myfile << row_ids;
       	 myfile.close();
    	} */	

	int idsend = 0;
	
 	node_id = nodeid_gen(0); 
	rows_node = row_ids.head(offset); 
	wheredata_node = wheredata.head(offset);


	/*ofstream myfile4 ("data/wheredata_node.dat");
        if (myfile4.is_open())
        {
         myfile4 << wheredata_node;
         myfile4.close();
        } */
	
   	 if (nodes > 1) {
		int count=0, count1=0; 
		for (int send = edison; send < world_size; send += edison) {
			//MPI_Send (&nodeid_gen(idsend), 1, MPI_INT, send, 121, MPI_COMM_WORLD);
        		MPI_Send (&row_ids(offset), chunksize, MPI_INT, send, 122, MPI_COMM_WORLD);
			MPI_Send (&wheredata(offset), chunksize, MPI_INT, send, 123, MPI_COMM_WORLD);
			idsend++;
			offset+=chunksize;
			
 		}	
	 }
  
    }

    if (nodes > 1) {
	//MPI_Recv (&node_id, 1, MPI_INT, 0, 121, MPI_COMM_WORLD, &status);
	MPI_Recv (wheredata_node.data(), offset, MPI_INT, 0, 123, MPI_COMM_WORLD, &status); 
	MPI_Recv (rows_node.data(), offset, MPI_INT, 0, 122, MPI_COMM_WORLD, &status);
	cout << "Received all 2 data" << endl;
    }

 
    /*if (world_rank == 24) {

	ofstream myfile4 ("data/2_wheredata_node.dat");
        if (myfile4.is_open())
        {
         myfile4 << wheredata_node;
         myfile4.close();
        }


	ofstream myfile5 ("data/2_rows_node.dat");
        if (myfile5.is_open())
        {
         myfile5 << rows_node;
         myfile5.close();
        }
	
    }*/
 
    MatrixXf d(m_node, n+1), d1(m_node, n), X_i(m_node, n+1); 
    VectorXf d2(m_node); 

    if(world_rank == MASTER)
                start_loadTime = MPI_Wtime();

    /* Create hyperslab selection and getting data into X,y from the hdf5 file */ 
	

    if (MPI_COMM_NULL != comm_ROOT) { 

	
	cout << "inside load area, loading.................." << endl;

	hid_t       file_id, dset_id, dset_id_y;         /* file and dataset identifiers */
	hid_t       filespace, memspace, filespace_y, memspace_y;      /* file and memory dataspace identifiers */
	hsize_t     count[2], count_y[1];                 /* hyperslab selection parameters */
    	hsize_t     offset[2], offset_y[1];
	hid_t       plist_id;                 /* property list identifier */
        herr_t      status;

	/* 
     	* Set up file access property list with parallel I/O access
     	*/

	cout << "before fapl" << endl; 

     	plist_id = H5Pcreate(H5P_FILE_ACCESS);
     	status =  H5Pset_fapl_mpio(plist_id, comm_ROOT, MPI_INFO_NULL);


	cout << "after fapl" << endl; 
    	/*
     	* Create a new file collectively and release property list identifier.
     	*/
    	file_id = H5Fopen(INPUTFILE, H5F_ACC_RDONLY, H5P_DEFAULT);

	/*
     	* Create the dataset with default properties and close filespace.
     	*/
    	dset_id = H5Dopen2(file_id, DATASET_X, H5P_DEFAULT);
	dset_id_y = H5Dopen2(file_id, DATASET_Y, H5P_DEFAULT);
	/* 
     	* Each process defines dataset in memory and writes it to the hyperslab
     	* in the file.
     	*/ 

    	count[0] = m_node;
    	count[1] = n;
    	offset[0] = root_rank * count[0];
    	offset[1] = 0;
    	memspace = H5Screate_simple(2, count, NULL);

	count_y[0] = m_node; 
	offset_y[0] = root_rank * count_y[0]; 
	memspace_y = H5Screate_simple(1, count_y, NULL);

	cout << "setting up memory stage" << endl; 

    	/*
     	* Select hyperslab in the file.
     	*/
    	filespace = H5Dget_space(dset_id);
    	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

	filespace_y = H5Dget_space (dset_id_y);
	H5Sselect_hyperslab(filespace_y, H5S_SELECT_SET, offset_y, NULL, count_y, NULL);

	cout << "Passed hyperslab stage" << endl; 

	/*
     	* Create property list for collective dataset write.
     	*/
   	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    	status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      H5P_DEFAULT, d1.data());

	cout << "loaded d1" << endl; 

	plist_id = H5Pcreate(H5P_DATASET_XFER);
        status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dread(dset_id_y, H5T_NATIVE_FLOAT, memspace_y, filespace_y,
                      H5P_DEFAULT, d2.data());


	cout << "Loaded d2" << endl; 

	if (root_rank == 0) { 
		end_loadTime = MPI_Wtime() - start_loadTime; 
		cout << "Load time: " << end_loadTime << endl;  
		start_commTime = MPI_Wtime(); 
	}

	d << d1,d2;   

	/* Communicate the data between the comm_ROOT by shuffled vector */ 

	/*for (int i=0; i < rows_node.size(); i++) {
	
		if (admm_rank == wheredata_node(i)) 
			X_i.row(i) = d.row(rows_node(i)); 

		else
			MPI_Sendrecv (d.row(rows_node(i)).data(), d.cols(), MPI_FLOAT, wheredata_node(i), 
				i, X_i.row(i).data(), X_i.cols(), MPI_FLOAT, node_id, i, MPI_COMM_WORLD, mpi_status);	
	
		}*/

	for (int i =0; i < rows_node.size(); i++) {
		if (root_rank == wheredata_node(i))
			 X_i.row(i) = d.row(rows_node(i));

		else
			MPI_Sendrecv (d.row(rows_node(i)).data(), d.cols(), MPI_FLOAT, wheredata_node(i),
                               	i, X_i.row(i).data(), X_i.cols(), MPI_FLOAT, root_rank, i, MPI_COMM_WORLD, mpi_status);

		}

	}

 	 

    if (world_rank == MASTER) {
	end_commTime = MPI_Wtime() - start_commTime;
	cout << "Communication time is: " << end_commTime << endl; 
		 /*ofstream myfile ("X.dat");
    		if (myfile.is_open())
    		{
        		myfile << X_i.row(1);
        		myfile.close();
    		} */
    } 


    if (world_rank == MASTER)
	start_distTime = MPI_Wtime(); 

    if (nodes > 1) 
    	MPI_Bcast (X_i.data(), m_node*(n+1), MPI_FLOAT, 0, comm_ADMM); 

    else
	 MPI_Bcast (X_i.data(), m_node*(n+1), MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (world_rank == MASTER) {
	end_distTime = MPI_Wtime() - start_distTime; 
	cout << "Data was distributed in " << end_distTime << endl; }

    VectorXf y; 
    MatrixXf X; 
    y  = X_i.rightCols(1); 
    X = X_i.topLeftCorner(m_node, n);  

    /*if (world_rank == 2) { 
    ofstream myfile9 ("data/X_20.dat");
    if (myfile9.is_open())
    {
        myfile9 << X;
        myfile9.close();
    }

   }*/ 

    if (world_rank == MASTER) {
	cout << "BoLBO analysis initialized" << endl; 
	cout << "------------------------------" << endl; 
	cout << "*No. of processes\t\t" << world_size << endl; 
	cout << "*No. of lambda elements\t\t" << nMP << endl; 
	cout << "*Max. no of iterations\t\t" << maxBoot << endl; 
	cout << "*No. of model dimensions\t" << n << endl; 
	cout << "*No. of model samples\t\t" << m << endl; 
	start_compTime = MPI_Wtime(); 
    }

   /*__________________________________
   MODEL SELECTION
   __________________________________*/


    VectorXd lamb0(nMP);
    //double my_lamb0; 

    
    lamb0.setLinSpaced(nMP, -3, 3);
    
    for(int i=0; i<nMP; i++)
    	lamb0(i) = pow(10,lamb0(i));
    
	
    //MPI_Scatter (lamb0.data(), 1, MPI_DOUBLE, &my_lamb0, 1, MPI_DOUBLE, 0, comm_ADMM); 

   
    double my_lamb0 = lamb0(mp_id); 

    VectorXf my_B0(n), R2m0(maxBoot*nMP), my_R2m0(1);
    MatrixXf B0(maxBoot*nMP,n), R2m0_m (maxBoot*nMP, 1); 
//    my_B0.setZero(); 
    R2m0.setZero(); my_R2m0.setZero(); B0.setZero(); R2m0_m.setZero();

    if (world_rank == MASTER)
	start_las1Time = MPI_Wtime();

    int m_train = round(rndfrctL*m_node); 

    
    VectorXi shuf_m, train_ids, ids;

    shuf_m.setLinSpaced(m_node, 0, m_node); 
    random_shuffle (shuf_m.data(), shuf_m.data()+shuf_m.size());
    train_ids = shuf_m.head(m_train); 

    ids = train_ids.segment(admm_rank, train_ids.size()-admm_rank);

    if (admm_rank == 0) {ids = shuf_m; }

    int my_m = ids.size(); ;

    VectorXi seeds(maxBoot);
    srand(seed); 
    for (int i=0; i<maxBoot; i++) {
 	int num = (rand() % 9999) + 1;
	seeds(i) = num; 
    }

    seed = seeds(bt_id); 
    srand(seed); 

    MatrixXf X_train;
    VectorXf y_train;


     if(admm_rank == 0) {
 	X_train = X.topRows(m_train); 
 	y_train = y.head(m_train); 	

      }

     else {
 	X_train = X.middleRows(admm_rank, m_train);
 	y_train = y.segment(admm_rank, m_train); 
     }

    VectorXf y_las; 
    y_las = y_train.array()-y_train.mean();
   
   /*ofstream myfile7 ("data/X_train.dat");
    if (myfile7.is_open())
    {
        myfile7 << X_train;
        myfile7.close();
    }*/
 
    //train
    my_B0 = lasso.lasso_admm (X_train, y_las, my_lamb0, comm_ADMM); 

//    if (world_rank == MASTER) { 
   // cout << "Printing my_B0" << endl; 
    /*ofstream myfile11 ("data/my_B0.dat");
    if (myfile11.is_open())
    {
        myfile11 << my_B0;
        myfile11.close();
    }*/
//}


    	//test
    MatrixXf X_test; 
    VectorXf y_test, yhat, y_mean; 
    //X_test = X.bottomRows(X.rows()-m_train); 
    //y_test = y.tail(y.size()-m_train); 
 
    /*ofstream myfile7 ("data/m_train.dat");
        if (myfile7.is_open())
        {
         myfile7 << m_train;
         myfile7.close();
        }
	
    ofstream myfile8 ("data/X_test.dat");
        if (myfile8.is_open())
        {
         myfile8 << X_test;
         myfile8.close();
        }*/

    //y_mean = y_test.array() - y_test.mean();
	
    //double sum  = y_test.sum();

   //cout << "\ny_mean\n" << y_test.colwise().mean() << endl; 
    /*ofstream myfile7 ("data/sum.dat");
        if (myfile7.is_open())
        {
         myfile7 << y_test.sum();
         myfile7.close();
        } */

    /*ofstream myfile ("data/y_test.dat");
    if (myfile.is_open())
    {
        myfile << y_test;
        myfile.close();
    }*/
 
//    if (admm_rank == 0) {
	X_test = X.bottomRows(X.rows()-m_train); 
    	y_test = y.tail(y.size()-m_train); 
    	yhat = X_test * my_B0; 
    	
	//cout << "\ny_mean\n" << y_test.mean() << endl;
	/*ofstream myfile1 ("data/yhat.dat");
    	if (myfile1.is_open())
    	{
       	 myfile1 << yhat;
       	 myfile1.close();
    	}*/

	y_mean = y_test.array() - y_test.mean();

	/*ofstream myfile1 ("data/y_mean.dat");
        if (myfile1.is_open())
        {
         myfile1 << y_mean;
         myfile1.close();
        }*/	

    	double r = num.pearson(yhat, y_mean);
    	my_R2m0 << r*r; 
	//cout << "\nyhat(10)\n" << yhat(10) << "r " << r << " rank = " << world_rank <<  endl; 
//    }

    //cout << "Before admm_rank == 0 line 524"  << endl;

    if (world_rank == MASTER) {
	end_las1Time = MPI_Wtime() - start_las1Time;
	start_bca1Time = MPI_Wtime(); 
    } 
  
    
    //MPI_Barrier (comm_ADMM);
    MPI_Gather (my_B0.data(), my_B0.size(), MPI_FLOAT, B0.data(), my_B0.size(), MPI_FLOAT, 0, comm_ADMM);
    MPI_Gather (my_R2m0.data(), 1, MPI_FLOAT, R2m0.data(), 1, MPI_FLOAT, 0, comm_ADMM);

    /*if (admm_rank == 0) {
	ofstream myfile1 ("data/R2m0.dat");
        if (myfile1.is_open())
        {
         myfile1 << R2m0;
         myfile1.close();
        }	

   }*/

     //cout << "Passed admm_rank == 0" << endl;
    MPI_Barrier (MPI_COMM_WORLD); 

    VectorXd lambL(nMP);

    if (world_rank == MASTER) {
	R2m0_m << R2m0; 
	R2m0_m.resize(maxBoot, nMP); 
	end_bca1Time = MPI_Wtime()-start_bca1Time;
	VectorXf Mt(nMP);
	VectorXi Lids;
	Mt = R2m0_m.colwise().mean() * 1e4; 
	transform (Mt.data(), Mt.data()+Mt.size(), Mt.data(), floorObject) ;
	Lids = num.where(Mt, Mt.maxCoeff()); 
	
	//cout << "R2m0\n" << R2m0 << endl;

	/*ofstream myfile1 ("data/R2m0_m.dat");
        if (myfile1.is_open())
        {
         myfile1 << R2m0_m;
         myfile1.close();
        }*/

	//cout << "Mt" << Mt << "\tMax Coeff" << Mt.maxCoeff() << endl;
	//cout << "LIDS\n" << Lids << endl;
	
	float v = lamb0(Lids(floor(Lids.size()/2)));
	float dv = pow(10, floor(log10(v)-1));

	//cout << "v= " << v << "\tdv = " << dv << endl;; 
	lambL.setLinSpaced(nMP, v-5*dv, v+5*dv);
    }

    MPI_Bcast(lambL.data(), nMP, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double my_lambL = lambL(mp_id);

//    if (world_rank==0)
//i    	cout << "lambL\n" << lambL << endl;

    VectorXf my_B(n), R2m(maxBoot*nMP), my_R2m(1);
    MatrixXf B(maxBoot*nMP,n), R2m_m(maxBoot*nMP,1);
    my_B.setZero(); R2m.setZero(); my_R2m.setZero(); B.setZero(); R2m_m.setZero();

    /*-------LASSO--------*/

    if (world_rank == MASTER) 
	start_las2Time = MPI_Wtime(); 

    //train
    my_B = lasso.lasso_admm (X_train, y_las, my_lambL, comm_ADMM);

//    if (admm_rank == 0) {
	X_test = X.bottomRows(X.rows()-m_train); 
    	y_test = y.tail(y.size()-m_train); 
    	//test
    	yhat = X_test * my_B; 
	/*ofstream myfile2 ("yhat.dat");
         if (myfile2.is_open())
         {
                 myfile2 << yhat;
                 myfile2.close();
        }

	ofstream myfile1 ("y_test.dat");
         if (myfile1.is_open())

         {
                 myfile1 << y_test.array()-y_test.mean();
                 myfile1.close();
        }*/

   	r = num.pearson(yhat, y_test.array()-y_test.mean());
    	my_R2m << r*r;
	//cout << "my_R2m= " << my_R2m << endl;
//    } 

 

    if (world_rank == MASTER) {
	end_las2Time = MPI_Wtime() - start_las2Time; 
	start_bca2Time = MPI_Wtime();
    }

    //cout << "in this loop" << endl; 
    //MPI_Barrier (comm_ROOT); 
    MPI_Gather (my_B.data(), n, MPI_FLOAT, B.data(), n, MPI_FLOAT, 0, comm_ADMM);
    MPI_Gather (my_R2m.data(), 1, MPI_FLOAT, R2m.data(), 1, MPI_FLOAT, 0, comm_ADMM);
   
   /*if (admm_rank == 0) {
        ofstream myfile1 ("data/R2m.dat");
        if (myfile1.is_open())
        {
         myfile1 << R2m;
         myfile1.close();
        }       

   } */ 

 

    MatrixXf sprt(nMP, n);
    //sprt.setZero();

    if (world_rank == MASTER) { 
	//R2m_m << R2m; 
	//R2m_m.resize(maxBoot, nMP);
	end_bca2Time = MPI_Wtime() - start_bca2Time;
	int b_row = 0;
	sprt.fill(NAN); 
	VectorXf intv, tmp_ids;

	for (int i=0; i < nMP; i++) {
		for (int j=0; j<nbootS; j++) { 
			tmp_ids = num.where_neq(B.row(b_row),0); 
	
			if (i == 0) { 

				
				/*ofstream myfile1 ("data/tmp_ids.dat");
        			 if (myfile1.is_open())
         			{
		                 myfile1 << tmp_ids;
                		 myfile1.close();
        			}*/

				intv = tmp_ids;

				/*ofstream myfile11 ("data/intv.dat");
                                 if (myfile11.is_open())
                                {
                                 myfile11 << intv;
                                 myfile11.close();
                                }*/				
			}

			intv = num.intersect1d(intv,tmp_ids);
			b_row++;		
		}
		sprt.row(i).head(intv.size()) = intv; 
	}

	/*ofstream myfile11 ("data/intv.dat");
                                 if (myfile11.is_open())
                                {
                                 myfile11 << intv;
                                 myfile11.close();
                                }*/


	/*ofstream myfile11 ("data/sprt.dat");
                                 if (myfile11.is_open())
                                {
                                 myfile11 << sprt;
                                 myfile11.close();
                                } */

	/*Sanity Check*/
	if (b_row != maxBoot*nMP) {cout << "\nCheck b_row iterator\n" << endl; }
    }
	
   //cout << "sprt (" << sprt.rows() << ":" << sprt.cols() << ")" << endl;

    /*ofstream myfile11 ("data/sprt.dat");
                                 if (myfile11.is_open())
                                {
                                 myfile11 << sprt;
                                 myfile11.close();
                                } */

    MPI_Bcast(sprt.data(), nMP * n, MPI_FLOAT, 0, comm_ADMM); 


    /*__________________________________
    MODEL ESTIMATION
    __________________________________*/

    m_frac = round (cvlfrct*m_node); 
    L_frac = round (rndfrct*m_frac); 


    MatrixXf Bgd(nrnd, n), rsd(nrnd, m_node-m_frac); 
    VectorXf R2(nrnd), bic(nrnd); 

    for (int cc=0; cc < nrnd; cc++) {
	VectorXf my_Bgols_B(n), my_Bgols_R2m(1); 
	my_Bgols_B.setZero(); my_Bgols_R2m.setZero(); 

	MatrixXf Bgols_B(maxBoot*nMP, n), Bgols_R2m_m(maxBoot*nMP, 1); 
	VectorXf Bgols_R2m(maxBoot*nMP);
	Bgols_B.setZero(); Bgols_R2m.setZero(); 

	VectorXi inds, L, T;
	inds.setLinSpaced(m_node, 0, m_node); 
	random_shuffle (inds.data(), inds.data()+inds.size()); 

	L = inds.head(m_frac); 
	//test set for final evaluation
	T = inds.tail(inds.size()-m_frac); 

	//train set is divided into training and testing set
	inds.setLinSpaced(m_frac, 0, m_frac); 
	random_shuffle(inds.data(), inds.data()+inds.size()); 
	VectorXi train, test; 
	train = inds.head(L_frac);
	test = inds.tail(inds.size()-L_frac); 

	VectorXi L_train;
	L_train = num.shuffles(L, train); 

	//TODO: if Cross Validation is across nodes do sendrecv
 
	if (admm_rank!=0) { train = num.select(L_train, ids);}
	else 		  { train = num.select_not(train_ids, L_train); }



	/****Select Support****/

	VectorXf sprt_row;  
	VectorXi sprt_ids, zdids, arange;

	sprt_ids = sprt.row(mp_id).unaryExpr(ptr_fun(is_not_NaN)).cast<int>();
	arange.setLinSpaced(n, 0, n); 
	zdids = num.setdiff1d(arange, sprt_ids); 


	/*--------Linear Regression---------*/

	VectorXf rgstrct; 
	 

	if (!sprt_ids.isZero()) {
		MatrixXf X_tr, y_tr; 
		X_tr = num.shuffles(X, train);
		y_tr = num.shuffles(y, train);

		rgstrct = lasso.lasso_admm (X_tr, y_tr.array()-y_tr.mean(), 0, comm_ADMM); 


		//apply support

		for (int i = 0; i < sprt_ids.size(); i++)
			if (sprt_ids(i) != 0)
				my_Bgols_B(sprt_ids(i)) = rgstrct(i); 
					 
 
		for (int i=0; i<zdids.size(); i++)
			my_Bgols_B(zdids(i)) = 0; 

		
		MatrixXf X_L; VectorXf y_L;
		X_L = num.shuffles(X, L_train);
		y_L = num.shuffles(y, L_train);

		yhat = X_L * my_Bgols_B; 
		double r = num.pearson(yhat, y_L.array()-y_L.mean()); 
		my_Bgols_R2m << r*r;

		//cout << "my_Bgols_R2m " << my_Bgols_R2m << endl;
	}

	if (world_rank == MASTER && cc == 0) {
		end_olsTime = MPI_Wtime() - start_olsTime;
		start_bca3Time = MPI_Wtime();
	}
	
//	MPI_Barrier (comm_ROOT);
	MPI_Gather (my_Bgols_B.data(), n, MPI_FLOAT, Bgols_B.data(), n, 
			MPI_FLOAT, 0, comm_ADMM);
 	MPI_Gather (my_Bgols_R2m.data(), 1, MPI_FLOAT, Bgols_R2m.data(), 1, 
				MPI_FLOAT, 0, comm_ADMM);	

	if (admm_rank == 0) {

		Bgols_R2m_m << Bgols_R2m; 
		Bgols_R2m_m.resize(maxBoot, nMP); 

		if (cc == 0)
			end_bca3Time = MPI_Wtime()-start_bca3Time; 

		/***Bagging******/
		VectorXf v;
		VectorXi  ids_r, ids_c;

		if (bgdOpt == 1) {
			v = Bgols_R2m_m.rowwise().maxCoeff();

			//cout << "v in line 830\n" << v << endl;

			/*ids_r = num.where_mv_row(Bgols_R2m_m, v);
			ids_c  = num.where_mv_col(Bgols_R2m_m, v); */
			MatrixXf btmp = MatrixXf::Zero(nbootE, n); 

			/*ofstream myfile4 ("data/Bgols_R2m_m.dat");
    			if (myfile4.is_open())
    			{   
        			myfile4 << Bgols_R2m_m;
       				myfile4.close();
    			}*/	

			/*ofstream myfile5 ("data/Bgols_B.dat");
                        if (myfile5.is_open())
                        {   
                                myfile5 << Bgols_B;
                                myfile5.close();
                        } */

			for (int kk =0; kk < nbootE; kk++) {
				VectorXi ids_0, ids_kk;
				/*ids_0 = num.where(ids_r, kk);
				ids_kk = num.where(ids_c,ids_0(kk));*/

				ids_kk = num.where(Bgols_R2m, v(kk)); 

				 /* ofstream myfile5 ("data/ids_0.dat");
                               	 if (myfile5.is_open())
                                	{
                                        myfile5 << ids_0;
                                        myfile5.close();
                                	}*/

				  /*ofstream myfile6 ("data/ids_c.dat");
                                if (myfile6.is_open())
                                {
                                        myfile6 << ids_c;
                                        myfile6.close();
                                }*/

				//cout << "\nids_kk size= " << ids_kk.size() << endl;

				/*ofstream myfile4 ("data/ids_kk.dat");
	                        if (myfile4.is_open())
        	                {   
                	                myfile4 << ids_kk;
                        	        myfile4.close();
                        	} */

				btmp.row(kk) = Bgols_B.row(ids_kk(floor(ids_kk.size()/2))); 	

				
			}
			//cout << "Passed this loop.." << endl; 
			Bgd.row(cc) = num.median(btmp); 
		}

		else {
			VectorXf mean_Bgols_R2m;
			MatrixXf Bgols_B_median(Bgols_B.rows(), ids.size());
			mean_Bgols_R2m = Bgols_R2m_m.colwise().mean(); 
			float vv = mean_Bgols_R2m.maxCoeff();
			ids_r = num.where(mean_Bgols_R2m, vv);

			for (int jj=0; jj<ids.size(); jj++)
				Bgols_B_median.col(jj) = Bgols_B.col(ids(jj));
			Bgd.row(cc) = num.median(Bgols_B_median); 	
		}
		MatrixXf X_T; 
		VectorXf y_T;

		X_T = num.shuffles(X, T);
		y_T = num.shuffles(y,T);
		
		//cout << "Before product" << endl;
		/*ofstream myfile6 ("data/Bgd_row.dat");
                                if (myfile6.is_open())
                                {
                                        myfile6 << Bgd.row(cc).transpose();
                                        myfile6.close();
                                }*/
		//cout << "X_T size are (" << X_T.rows() << ":" << X_T.cols() << ")" << endl; 
		//cout << "Bgd row size " << Bgd.row(cc).size() << endl;  	

		yhat = X_T * Bgd.row(cc).transpose(); 
		double r = num.pearson(yhat, y_T.array()-y_T.mean());
		R2(cc)  = r*r;

		/*ofstream myfile6 ("data/y_hat.dat");
                                if (myfile6.is_open())
                                {
                                        myfile6 << yhat;
                                        myfile6.close();
                                }*/ 
		rsd.row(cc) = y_T.array()-y_T.mean()-yhat.array(); 
		bic(cc) = (m_node-m_frac) * log(rsd.row(cc).squaredNorm()/(m_node-m_frac)) + log(m_node-m_frac) * n;

		//cout << "Passed bic" << endl;  
	}
    }

   /* if (admm_rank == 0) {
        ofstream myfile1 ("data/R2.dat");
        if (myfile1.is_open())
        {
         myfile1 << R2;
         myfile1.close();
        }       

   }*/

    MPI_Barrier(MPI_COMM_WORLD); 
   
    if ( world_rank == MASTER)
     {

	//cout << "Inside printing......" << endl; 
    	end_compTime = MPI_Wtime() - start_compTime; 
    	end_ccloopTime = MPI_Wtime() - start_ccloopTime; 
	start_saveTime = MPI_Wtime();
	
    
	hid_t       outfile, dataset, dataspace, group_id, attribute_id;
    	herr_t      status; 

    	double data1[1], data2[1], data3[1], data4[1], data5[1], data6[1], data7[1], data8[1], data9[1], data10[1], data11[1], data12[1], data13[1], data14[1], data15[1], data16[1], data17[1], data18[1], data19[1], data20[1], data21[1], data22[1];
 
    	data1[0] = bgdOpt; 
    	data2[0] = nrnd; 
    	data3[0] = cvlfrct; 
    	data4[0] = rndfrct; 
    	data5[0] = rndfrctL; 
    	data6[0] = nbootE;
    	data7[0] = nbootS; 
    	data8[0] = nMP; 
    	data9[0] = seed; 
    	data10[0] = m;
    	data11[0] = n;
 
    	data12[0] = end_sizeTime; 
    	data13[0] = end_commTime; 
    	data14[0] = end_loadTime;  
    	data15[0] = end_distTime; 
    	data16[0] = end_compTime; 
    	data17[0] = end_las1Time; 
    	data18[0] = end_bca1Time; 
    	data19[0] = end_las2Time; 
    	data20[0] = end_bca2Time; 
    	data21[0] = end_olsTime; 
    	data22[0] = end_bca3Time; 
    
    
   	hsize_t  dimsfs[1], dimsf[2], dimst[3]; 

    	dimsfs[0] = 1;

    	outfile = H5Fcreate ("stats.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    	dataspace = H5Screate_simple (1, dimsfs, NULL);
    	dataset = H5Dcreate (outfile, "Units", H5T_IEEE_F32LE , dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	attribute_id = H5Acreate2 (dataset, "bdgOpt", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);

    	status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data1);

    	H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "nrnd", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data2);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "cvlfrct", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data3);

        H5Aclose(attribute_id);
 
        attribute_id = H5Acreate2 (dataset, "rndfrct", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data4);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "rndfrctL", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data5);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "nbootE", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data6);

        H5Aclose(attribute_id);
  
        attribute_id = H5Acreate2 (dataset, "nbootS", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data7);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "nMP", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data8);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "seed", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data9);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "m", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data10);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "n", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data11);

        H5Aclose(attribute_id);
        H5Dclose(dataset);
        H5Sclose(dataspace);

    	dataspace = H5Screate_simple (1, dimsfs, NULL);
        dataset = H5Dcreate (outfile, "Times", H5T_IEEE_F32LE , dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        attribute_id = H5Acreate2 (dataset, "sizeTime", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data12);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "commTime", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data13);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "loadTime", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data14);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "distTime", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data15);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "compTime", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data16);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "las1Time", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data17);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "bca1Time", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data18);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "las2Time", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data19);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "bca2Time", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data20);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "olsTime", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data21);

        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (dataset, "bca3Time", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, data22);

        H5Aclose(attribute_id);
        H5Dclose(dataset);
        H5Sclose(dataspace);

	group_id = H5Gcreate2(outfile, "/lasso", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	dimsfs[0] = nMP; 
    dataspace = H5Screate_simple (1, dimsfs, NULL);
        dataset = H5Dcreate (outfile, "lamb0",  H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lamb0.data());

    H5Dclose(dataset);
        H5Sclose(dataspace);

    dimsfs[0] = nMP; 
        dataspace = H5Screate_simple (1, dimsfs, NULL);
        dataset = H5Dcreate (outfile, "lambL",  H5T_IEEE_F64LE , dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lambL.data());

   	H5Dclose(dataset);
        H5Sclose(dataspace);
    	H5Gclose (group_id);

	group_id = H5Gcreate2(outfile, "/bolbo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	dimsf[0] = nMP;
    dimsf[1] = n;  
        dataspace = H5Screate_simple (2, dimsf, NULL);
        dataset = H5Dcreate (outfile, "sprt",    H5T_IEEE_F64LE , dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sprt.data());
	
	H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Gclose (group_id);
	H5Fclose(outfile);
    }

    if (admm_rank == 0) {
	
	hid_t       file_id, dset_id;         /* file and dataset identifiers */
    	hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    	hsize_t     dimsf[2];                 /* dataset dimensions */
        hsize_t     count[2];                 /* hyperslab selection parameters */
    	hsize_t     offset[2];
    	hid_t       plist_id;                 /* property list identifier */
    	herr_t      status;
	hsize_t     dimst[3], dimsfs[1];
	hsize_t     offset_t[3], offsets[1];
	hsize_t     count_t[3], counts[1];

	/* 
     	* Set up file access property list with parallel I/O access
     	*/
     	plist_id = H5Pcreate(H5P_FILE_ACCESS);
     	H5Pset_fapl_mpio(plist_id, comm_ROOT, MPI_INFO_NULL);

    	/*
     	* Create a new file collectively and release property list identifier.
     	*/
    	file_id = H5Fcreate(OUTPUTFILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    	H5Pclose(plist_id);


    	/*
     	* Create the dataspace for the dataset.
     	*/
    	dimsf[0] = nrnd * nodes;
    	dimsf[1] = n;
    	filespace = H5Screate_simple(2, dimsf, NULL);

    	/*
     	* Create the dataset with default properties and close filespace.
    	*/
    	dset_id = H5Dcreate(file_id, "Bgd", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    	H5Sclose(filespace);

    	/* 
     	* Each process defines dataset in memory and writes it to the hyperslab
     	* in the file.
     	*/
    	count[0] = dimsf[0]/nodes;
    	count[1] = dimsf[1];
    	offset[0] = node_id * count[0];
    	offset[1] = 0;
    	memspace = H5Screate_simple(2, count, NULL);

   	/*
     	* Select hyperslab in the file.
     	*/
    	filespace = H5Dget_space(dset_id);
    	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

	/*
     	* Create property list for collective dataset write.
     	*/
    	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    	status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, Bgd.data());


	dimsf[0] = maxBoot * nodes;
	dimsf[1] = nMP; 
        filespace = H5Screate_simple(2, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "R2m0", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /* 
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        count[0] = dimsf[0]/nodes;
        count[1] = dimsf[1];
        offset[0] = node_id * count[0];
        offset[1] = 0;
        memspace = H5Screate_simple(2, count, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, R2m0.data());


	dimsf[0] = maxBoot * nodes;
        dimsf[1] = nMP;
        filespace = H5Screate_simple(2, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "R2m", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /* 
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        count[0] = dimsf[0]/nodes;
        count[1] = dimsf[1];
        offset[0] = node_id * count[0];
        offset[1] = 0;
        memspace = H5Screate_simple(2, count, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, R2m.data());

	dimsf[0] = nrnd * nodes;
        dimsf[1] = m_node - m_frac;
        filespace = H5Screate_simple(2, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "rsd", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /* 
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        count[0] = dimsf[0]/nodes;
        count[1] = dimsf[1];
        offset[0] = node_id * count[0];
        offset[1] = 0;
        memspace = H5Screate_simple(2, count, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, rsd.data());


	dimsf[0] = nrnd * nodes;
        dimsf[1] = m_node - m_frac;
        filespace = H5Screate_simple(2, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "bic", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /* 
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        count[0] = dimsf[0]/nodes;
        count[1] = dimsf[1];
        offset[0] = node_id * count[0];
        offset[1] = 0;
        memspace = H5Screate_simple(2, count, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, bic.data());


	dimst[0] = maxBoot * nodes;
        dimsf[1] = nMP;
	dimst[2] = n;
        filespace = H5Screate_simple(3, dimst, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "B0", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /* 
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        count_t[0] = dimst[0]/nodes;
        count_t[1] = dimst[1];
	count_t[2] = dimst[2]; 
        offset_t[0] = node_id * count[0];
        offset_t[1] = 0;
	offset_t[2] = 0;
        memspace = H5Screate_simple(3, count_t, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, B0.data());


	dimst[0] = maxBoot * nodes;
        dimsf[1] = nMP;
        dimst[2] = n;
        filespace = H5Screate_simple(3, dimst, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "B", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /*
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        count_t[0] = dimst[0]/nodes;
        count_t[1] = dimst[1];
        count_t[2] = dimst[2]; 
        offset_t[0] = node_id * count[0];
        offset_t[1] = 0;
        offset_t[2] = 0;
        memspace = H5Screate_simple(3, count_t, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, B.data());


	dimsfs[0] = nrnd * nodes;
        filespace = H5Screate_simple(1, dimsfs, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "R2", H5T_NATIVE_FLOAT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /*
        * Each process defines dataset in memory and writes it to the hyperslab
        * in the file.
        */
        counts[0] = dimsfs[0]/nodes;
        offsets[0] = node_id * count[0];
        memspace = H5Screate_simple(1, counts, NULL);

        /*
        * Select hyperslab in the file.
        */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
        * Create property list for collective dataset write.
        */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, R2.data());

    }

    if (world_rank == MASTER) {
	end_saveTime = MPI_Wtime() - start_saveTime;         
	cout << "\nBolbo analysis completed" << endl; 
   	cout << "---------------------------------------" << endl;
   	cout << "Results stored in " << OUTPUTFILE << endl; 
   	cout << "Bolbo Times: " << endl;
	cout << "-size time: " << end_sizeTime << endl;
	cout << "-comm time: " << end_commTime << endl;
	cout << "-load time: " << end_loadTime << endl;
	cout << "-dist time: " << end_distTime << endl;
	cout << "-comp time: " << end_compTime << endl;
	cout << "\t-las1 time: " << end_las1Time << endl;
	cout << "\t-bca1 time: " << end_bca1Time << endl;
	cout << "\t-las2 time: " << end_las2Time << endl;
	cout << "\t-bca2 time: " << end_bca2Time << endl;
	cout << "\t-ols time: " << end_olsTime << endl;
	cout << "\t-bca3 time: " << end_bca3Time << endl;
	cout << "-save time: " << end_saveTime << endl;
    } 


    MPI_Finalize();

}
