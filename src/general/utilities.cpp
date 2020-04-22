#include "utilities.h"
#include "globheads.h"
#include "protos.h"

using namespace std;
typedef map<int, set<int> > Graphmap;

// Fills v with random values
void randomvec (vector<DataT> &v, int n) {
  /* fills v with random values */
  DataT x;
  int k;
  for (k=0; k<n; k++) {
    x = rand()%100;
    v[k] = x/100;
  }
}
//-----------------------------------------------------------------//


//Create a random sparse matrix
void random_sparse_mat(Mat &mat, int N, float sparsity){
    mat.row_size = N;
    mat.vec = vector<DataT>(N*N,0.);
    randomvec(mat.vec,(int) N*N*sparsity);
    std::random_shuffle(mat.vec.begin(),mat.vec.end());
}
//-----------------------------------------------------------------//



void graph_print(Graphmap& gmap) {
	for (auto const& node : gmap) {
		cout << "NODE:" << node.first << " CHILDREN:";
		for (auto child : node.second) {
			cout << child << " ";
		}
		cout << endl;
	}

}
//-----------------------------------------------------------------//



//-------------------------------------------------------------------------------------------------------------
//GRAPH UTILITIES
//
//
//
//TODO RMAT reader/generator (maybe in Python)


void read_snap_format(Graphmap& gmap, string filename)
{
/*		Read from edgelist to a graphmap
 *		TODO Error handling
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * filename: the file were the edgelist is stored
 *----------------------------------------------------------------------------
 * on return:
 * ==========
 * gmap   = the edgelist now into Graphmap format: map<int, set<int>>
 *----------------------------------------------------------------------------
 */

	gmap.clear();

	ifstream infile;

	//TODO handle reading error
	infile.open(filename);
	string temp = "-1";
	int current_node = -1, child;
	set<int> emptyset;

    // Ignore comments headers
    while (infile.peek() == '%') infile.ignore(2048, '\n');
    
	//TODO: check import success
	//read the source node (row) of a new line
	while (getline(infile, temp, ' ')) {
		current_node = stoi(temp);

		if (gmap.count(current_node)==0) { //new source node encountered
			gmap[current_node] = emptyset;
		}

		//get target node (column)
		getline(infile, temp);
		child = stoi(temp);
		gmap[current_node].insert(child);
	}
}

//check if a graphmap is well-numbered, complete and (optional) symmetric
int isProper(const Graphmap& gmap, bool mirroring) {
	//returns:
	//0 : everything ok
	//1 : there is a jump in node numeration
	//2 : inexistent children
	//3 : asymmetric link  (only if mirroring is on)

	int check = 0;
	int last_node = -1;

	for (auto const& x : gmap)
	{
		//check if there is a jump
		if (x.first > last_node + 1) { return 1; }
		last_node = x.first;

		//check node children
		for (auto const& elem : x.second)
		{
			//check if the children exists in the graph
			auto child = gmap.find(elem); //search result
			if (child == gmap.end()) {
				//inexistent children
				return 2;
			}
			//check for asymmetric links (only if mirroring is on)
			else if (mirroring) {
				if ((child->second).find(x.first) == (child->second).end()) {
					return 3;
				}
			}
		}
	}
	return 0;
}


//Make a graph undirected
void MakeUndirected(Graphmap& gmap) {

	for (auto x : gmap)
	{
		//check node children
		for (auto child : x.second) {
			(gmap[child]).insert(x.first);
		}
	}
}

//rename graph nodes so that they are consecutive integers 
//Quite costly. Probably worth checking with isProper first
void MakeProper(Graphmap& gmap) {

	map<int, int> new_name;
	int count = -1;
	int name = -1;
	set<int> tempset;

	//find correct names
	for (auto parent : gmap)
	{
		count++;
		name = parent.first;
		new_name[name] = count;
		cout << "name: " << name << " count: " << count << endl;
	}

	//change target names
	for (auto parent : gmap)
	{

		tempset.clear();
		for (auto child : parent.second) {
			tempset.insert(new_name[child]);
			cout << "child: " << child << " newname: " << new_name[child] << endl;
		}
		gmap[parent.first] = tempset;
	}

	//change source names
	for (auto n : new_name)
	{
		if (n.first != n.second) {
			gmap[n.second] = gmap[n.first];
			gmap.erase(n.first);
		}
	}
}

//TODO (not urgent)
//Conversion to edgelist and from CSR;
//weighted graphs
//options for unweighted graphs (e.g. only structure when calculating permutations)
//efficient implementation of undirected graphs (use symmetry where possible)

//Analysis of multiple graphs

//export as edgelist
void write_snap_format(Graphmap& gmap, string filename) {
	ofstream outfile;
	outfile.open(filename);
	int name;

	for (auto parent : gmap) {
		name = parent.first;

		for (auto child : parent.second) {
			outfile << name << " " << child << endl;
		}

	}
	outfile.close();
}

/*-------------------------------------------------------------------------------------------------------------
END OF GRAPH UTILITIES
*/

//Convert a sparse matrix into a Mat one
void convert_from_CSR(SparMat &spmt, Mat &mat){
    int n = spmt.n;
    mat.row_size = n;
    int idx;
    DataT val;
    mat.vec = vector<DataT>(n*n,0.);//initialize Mat's vec to 0
    
    //loop through rows
    for (int row = 0; row < n; row++){
        
        //loop through non-zeros in that row
        for (int i = 0; i< spmt.nzcount[row];i++){
            idx = row*n + spmt.ja[row][i];//position of nz value
            val = spmt.ma[row][i]; //nz value
            mat.vec[idx] = val; //put it in the proper position
        }
    }
    
}

//fill a CSR from three vectors
void fill_CSR(SparMat &spmt, int n, const vector<int> &vec_nzcount, const vector<int> &vec_ja,  const vector<DataT> &vec_ma){
    
    //declare the arrays to be put in the SparMat
    spmt.n = n;
    spmt.nzcount = new int[n];
    spmt.ja = new int*[n];
    spmt.ma = new DataT*[n];
    
    //populate the SparMat arrays using the vectors;
    std::copy(vec_nzcount.begin(),vec_nzcount.end(),spmt.nzcount);

    int nzs = 0; //total non-zeros written
    
    //loop through rows
    for (int row = 0; row < n; row++){
        //allocate ja,ma for that row
        spmt.ja[row] = new int[spmt.nzcount[row]];
        spmt.ma[row] = new DataT[spmt.nzcount[row]];
        
        //copy from vectors
        std::copy(vec_ma.begin() + nzs,vec_ma.begin() + nzs + spmt.nzcount[row],spmt.ma[row]);
        std::copy(vec_ja.begin() + nzs,vec_ja.begin() + nzs +  spmt.nzcount[row],spmt.ja[row]);
        nzs += spmt.nzcount[row];
    }
}


//Convert a Mat into a sparse matrix;
void convert_to_CSR(const Mat &mat, SparMat &spmt){
    vector<int> vec_ja, vec_nzcount;
    vector<DataT> vec_ma;
    int n = mat.row_size;
    setupCS(&spmt,n,1);//allocate memory
    
    //build the appropriate vectors from Mat;
    int i = 0;
    int row_count = 0;
    DataT temp;
    
    //loop through all elements in the Mat
    while ( i < pow(n,2) ){
        temp = mat.vec[i];
        
        //select non-zeros
        if (temp != 0.){
            vec_ma.push_back(temp);//value
            vec_ja.push_back(i % n);//column
            row_count++;
        }
        //check end of row
        if ((i+1) % n == 0){
            vec_nzcount.push_back(row_count);//nz in that row
            row_count = 0;
        }
        i++;
    }
    fill_CSR(spmt, n, vec_nzcount, vec_ja, vec_ma);
}

//Graphmap To CSR
//TODO TEST
void convert_to_CSR(const Graphmap &gmap, SparMat &spmt) {

	vector<int> vec_ja, vec_nzcount;
	vector<DataT> vec_ma;
	int n = gmap.size();
	setupCS(&spmt, n, 1);//allocate memory

	//build the appropriate vectors from Mat;
	for (auto node : gmap) {

		set<int> tempset = node.second; //adiacency list for the node.

		
		vec_ja.insert(vec_ja.end(), tempset.begin(), tempset.end()); //columns of nonzero elements = names of node children

		vector<DataT> tempvec(tempset.size(),1.); 
		vec_ma.insert(vec_ma.end(), tempvec.begin(), tempvec.end()); //entries = 1s. Graphmap are unweighted for now.
		
		vec_nzcount.push_back(tempset.size());//nonzero per row = number of node children
	}
	
	//use vectors to fill CSR
	fill_CSR(spmt, n, vec_nzcount, vec_ja, vec_ma);
}

void permute(SparMat &spmt, int* perm){
    //permutes a matrix in CSR form
    int n = spmt.n;
    if (permute(spmt.nzcount,perm,n) or permute(spmt.ja,perm,n) or permute(spmt.ma,perm,n)) cout << "ATTENTION: COULD NOT PERFORM PERMUTATION" << endl;
}

void matprint(const SparMat &spmt){
    cout << "\n  ********************************** \n PRINTING THE CSR MATRIX" << endl;
    int n = spmt.n;
    int sparse_count = spmt.nzcount[n-1];
    
    cout << "n = " << n << endl;

    cout << endl << "printing nzcount" << endl;
    std::copy(spmt.nzcount, spmt.nzcount + n, std::ostream_iterator<int>(std::cout, " "));

    cout << endl << "printing ja" << endl;
    for (int row = 0; row < n; row++){
        for (int i = 0; i <spmt.nzcount[row]; i++){
            cout << spmt.ja[row][i] << " ";
        }
        cout << endl;
    }
    
    cout << endl << "printing ma" << endl;
    for (int row = 0; row < n; row++){
	cout << setprecision(5);
        for (int i = 0; i <spmt.nzcount[row]; i++){
            cout << spmt.ma[row][i] << " ";
        }
        cout << endl;
    }
    
    cout << "done" << endl;
}

void matprint(const Mat &mat){
    cout << "PRINTING THE MATRIX" << endl;
    int N = mat.row_size;
    for (int i = 0; i < N*N; i++){
        if (i%N == 0){
            cout << endl;
	    cout << setprecision(5);
        }
        cout << mat.vec[i] << " ";
    }
    cout <<endl<< "finished printing the matrix"<<endl;

}

void matprint(const DataT *mat, const int rows, const int cols, const bool transpose)
{

        for(int i = 0; i<rows;i++){
                for (int j = 0; j < cols; j++){
                        if (transpose){
                                std::cout<< mat[j*rows + i]<< " ";
                        }
                        else{
                                std::cout<< mat[i*cols + j]<< " ";
                        }
                }
                std::cout << std::endl;
        }

}


void matprint(const VBSparMat &vbmat){
    int N = vbmat.n, *bsz = vbmat.bsz;
    int Lsz,Hsz,col;
    cout << "PRINTING a VB MATRIX" << endl;
    cout<<"DIVIDED IN ="<< N <<endl;
    
    for(int i = 0; i < N; i++ ) {
    	cout<<"bsz["<<i+1<<"] = "<<bsz[i+1]<<endl;
        Hsz = bsz[i+1] - bsz[i];
        cout<<"nzcount["<<i<<"] = "<<vbmat.nzcount[i]<<endl;
        
	for(int j = 0; j<vbmat.nzcount[i]; j++){
            col = vbmat.ja[i][j];
            cout<<"ja["<<i<<"]["<<j<<"] = "<<col<<endl;
            Lsz = bsz[col+1] - bsz[col];
            cout<<"BLOCK size="<<Hsz<<"x"<<Lsz<<endl;
            
	    for(int k_row = 0; k_row < Hsz; k_row++){
		for( int k_col = 0; k_col < Lsz; k_col++){
                	cout<<setprecision(5)<<vbmat.ba[i][j][Hsz*k_col + k_row]<<" ";
		}
		cout <<endl;

            }
            cout<<endl;
        }
    }
}


//read an unweighted graph into a SparMat from a edgelist file
//edges must be ordered by source
void read_snap_format(SparMat &spmt,
	string infilename) {

	ifstream infile;
	infile.open(infilename);
	string temp = "-1";

	int row = -1, col;
	int row_count = 1;
	int sparse_count = 0;
	int tot_count = 0;
	vector<int> vec_ja, vec_nzcount;
	vector<DataT> vec_ma;
	int last_node = -1;

	//read the source node (row) of a new line
	while (getline(infile, temp, ' ')) {

		//count total entries in the matrix
		sparse_count++;

		if (stoi(temp) != last_node) { //conta i nodi diversi

			//finalize last row if there is one
			if (last_node != -1) {
				vec_nzcount.push_back(row_count);
				row_count = 1;
			}

			//insert blank rows for disconnected nodes
			for (int i = stoi(temp) - 1; i > last_node; i--) {
				vec_nzcount.push_back(0);
				tot_count++;
			}

			tot_count++;
			last_node = stoi(temp);
		}
		else{
			row_count++;
		}

		//get target node (column)
		getline(infile, temp);
		col = stoi(temp);
		vec_ma.push_back(1.);
		vec_ja.push_back(col);
	}

	//finalize last row
	vec_nzcount.push_back(row_count);
    
    fill_CSR(spmt, tot_count, vec_nzcount, vec_ja, vec_ma);
}

//read a Matrix Market sparse matrix from file into a SparMat
//extremely inefficient, first build a Mat and then converts to CSR
void read_mtx_format(SparMat &spmt, string infilename){
    ifstream file(infilename);
    int num_row, num_col, num_lines;
    
    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;
    
    Mat mat;
    mat.row_size = num_row;
    vector<DataT> temp_mat(num_row*num_row, 0.0);
    
    // fill the matrix with data
    for (int l = 0; l < num_lines; l++)
    {
        DataT data;
        int row, col;
        file >> row >> col >> data;
        temp_mat[(row -1) + (col -1) * num_row] = data;
    }
    
    mat.vec = temp_mat;
    convert_to_CSR(mat,spmt);
    
}



/*extracts some features of a Block Matrix vbmat:
Msize: side length of the matrix
Bnum: number of non-zero blocks
BLsize,BHsize: height and length of each block
Bsparsity: sparsity of each block
*/
void extract_features(const vbsptr vbmat, int& Msize, int &Bnum, vector<int>& BLsize, vector<int>& BHsize, vector<int>& Bsparsity){
    int N = vbmat->n, *bsz = vbmat->bsz;
    int Lsz,Hsz,col;
    double tempSparse;
    Msize = bsz[N];
    Bnum = 0;
    //loop vertically through block rows
    for(int i = 0; i < N; i++ ) {
        Hsz = bsz[i+1] - bsz[i];
        Bnum += vbmat->nzcount[i];
        //loop horizontaly through block columns
        for(int j = 0; j<vbmat->nzcount[i]; j++){
            col = vbmat->ja[i][j];
            Lsz = bsz[col+1] - bsz[col];
            BLsize.push_back(Lsz);
            BHsize.push_back(Hsz);
			tempSparse = 0;
            //loop (column-wise) through elements in the block
            for(int k = 0; k < Lsz*Hsz; k++){
                if (vbmat->ba[i][j][k] != 0.){
                    tempSparse += 1;
                }
            }
            Bsparsity.push_back(tempSparse);
        }
    }
}

//TODO output to string, not CSV
void features_to_CSV(vbsptr vbmat, ofstream& CSV_out, int verbose = 0){
    int Msize,Bnum;
    int BlockTotal = 0;
	int SparseTotal = 0;
    vector<int>BLsize,BHsize,Bsparsity;
    cout<<"EXTRACTING AND PRINTING FEATURES"<<endl<<endl;
    extract_features(vbmat, Msize , Bnum , BLsize, BHsize, Bsparsity);
    
    for(int i = 0; i != BLsize.size(); ++i){
            BlockTotal += BLsize[i]*BHsize[i];
            SparseTotal += Bsparsity[i];
        }

	if (verbose == 1) {
		cout << "Matrix size = " << Msize << "X" << Msize << endl;
		cout << "Original Sparsity = " << (double)SparseTotal / (Msize*Msize) << endl;
		cout << "Divided in = " << vbmat->n << endl;
		cout << "Number of nonzero blocks = " << Bnum << " out of " << vbmat->n*(vbmat->n) << endl;
		cout << "Average height of nonzero blocks = " << mean(BHsize) << " +- " << stdv(BHsize)<< endl;
		cout << "Average length of nonzero blocks = " << mean(BLsize) << " +- " << stdv(BLsize)<< endl;
		cout << "New Area fraction = " << (double)BlockTotal / (Msize*Msize) << endl;
		cout << "Average nonzeros in blocks = " << mean(Bsparsity) << " +- " << stdv(Bsparsity) << endl;
		cout << "New Sparsity = " << (double)SparseTotal / BlockTotal << endl;
	}

	//print on CSV

	CSV_out << Msize << "," << (double)SparseTotal / (Msize*Msize) << ","
		<< vbmat->n << "," << Bnum << "," 
		<< mean(BHsize) << "," << stdv(BHsize) << ","
		<<  mean(BLsize)<< "," << stdv(BLsize) << ","
		<< (double)BlockTotal / (Msize*Msize) << ","
		<< mean(Bsparsity) << "," << stdv(Bsparsity) << ","
		<< (double)SparseTotal / BlockTotal << endl;
	
}

//reorder a CSR matrix spmt and generate a Block Sparse Matrix
//TODO return permutation perm for further use
int make_sparse_blocks(SparMat &spmat, VBSparMat &vbmat, double eps){
    int nBlock;
    int *nB = NULL, *perm = NULL;
    double *t_hash = NULL, *t_angle = NULL;
    
    if (init_blocks(&spmat, &nBlock, &nB, &perm, eps, t_hash, t_angle) != 0) {
        cout << "ERROR: COULD NOT CREATE PERMUTATION. Try with another epsilon" << endl;
        return -1;
    }


    
    permute(spmat,perm);    
   


    int ierr = csrvbsrC(1, nBlock, nB, &spmat, &vbmat); //convert SparMat spmat into VBSparMat vbmat, using block structure found by init_block
    if (ierr != 0){
        cout << "error occurred while creating block matrix" << endl;
        return -1;
    }

}


int random_sparse_blocks_mat(Mat &mat, int N, int n_block, float block_k, float k){


	//initialize Mat	
	mat.row_size = N;
    	mat.vec = vector<DataT>(N*N,0.);
    	int row_count = 0;
	
	//find number of nonzero entries in matrix and blocks
	int nzcount = (int) k*N;
	float k_inside_blocks = k/block_k;
	if (k_inside_blocks >1) {
		block_k = 1./k;
		k_inside_blocks = 1.;
		cout << "WARNING: block sparsity must be greater than 1/k. Changing block_k to "<< block_k <<endl;
	}	
	
	//fix nonzero blocks
	int nzblocks = (int) n_block * n_block * block_k; //nonzero block number
	vector<int> blocks = vector<int> ( n_block*n_block, 0); //one element per block, 0 if empty, 1 otherwise
	std::fill (blocks.begin(),blocks.begin()+nzblocks,1); 
	std::random_shuffle(blocks.begin(),blocks.end());
	
	int block_dim = (int) N/n_block;
	if (block_dim == 0) {
		block_dim = 1;
		n_block = N;
		cout << "WARNING: block number must be lower or equal to N. Changing n_block to N"<< n_block <<endl;
	}
	int nz_in_block = (int) block_dim * block_dim * k_inside_blocks; //nonzeros in a block
	vector<DataT> block_vals;
		
	//put nonzerovalues in the Mat
	for(int ib = 0; ib < n_block; ib++){//iterate through block rows
		for (int jb = 0; jb < n_block; jb++){ //iterate through block columns
			if(blocks[ib*n_block + jb] != 0){
				//if block is nonempty, put random values in it;
				block_vals = vector<DataT>(block_dim * block_dim, 0);
        			randomvec(block_vals,(int) nz_in_block);
				std::random_shuffle(block_vals.begin(),block_vals.end());
				//put the element in the vector representing the block in the right positions in the Mat
				for (int i = 0; i < block_dim; i++){
					for (int j = 0; j <block_dim; j++){
						int row = ib*block_dim + i;
						int col = jb*block_dim + j;
						mat.vec[row*N + col] = block_vals[i*block_dim + j];
					}
				}
			}
		}
	}


}

//print array
template <class myType>
void arrprint(myType *arr, int len){
    for(int i = 0; i < len; i++){
        cout<<arr[i]<<" ";
    }
    cout<<endl;
}

void convert_to_col_major(DataT *X, DataT *Y, const int rows, const int cols){
	for (int i=0;i<rows;i++){
		for (int j=0; j<cols; j++){
			Y[rows*j + i] = X[cols*i + j];
		}
	}

}

void convert_to_row_major(DataT *X, DataT *Y, const int rows, const int cols){
        for (int j=0;j<cols;j++){
                for (int i=0; i<rows; i++){
                        Y[cols*i + j] = X[rows*j + i];
                }
        }

}

bool are_equal(const DataT *X,const DataT* Y,const int n, const double eps){
	for (int i = 0; i < n; i++){
		if(abs(X[i] - Y[i]) > eps) {
//			cout<<"pos "<< i<< " : " << X[i]<< " !=  " << Y[i] <<endl;
			return false;
		}
	}
	return true;
}
