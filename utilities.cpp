#include "utilities.h"
#include "globheads.h"
#include "protos.h"

using namespace std;
typedef map<int, set<int> > Graphmap;

// Fills v with random values
void randomvec (vector<double> &v, int n) {
  /* fills v with random values */
  double x;
  int k, seed = 4321;
  srand(seed);
  for (k=0; k<n; k++) {
    x = rand()%100;
    v[k] = x/100;
  }
}
//-----------------------------------------------------------------//


//Create a random sparse matrix
void random_sparse_mat(Mat &mat, int N, double sparsity){
    mat.row_size = N;
    mat.vec = vector<double>(N*N,0.);
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

//TODO
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
    double val;
    mat.vec = vector<double>(n*n,0.);//initialize Mat's vec to 0
    
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
void fill_CSR(SparMat &spmt, int n, const vector<int> &vec_nzcount, const vector<int> &vec_ja,  const vector<double> &vec_ma){
    
    //declare the arrays to be put in the SparMat
    spmt.n = n;
    spmt.nzcount = new int[n];
    spmt.ja = new int*[n];
    spmt.ma = new double*[n];
    
    //populate the SparMat arrays using the vectors;
    std::copy(vec_nzcount.begin(),vec_nzcount.end(),spmt.nzcount);

    int nzs = 0; //total non-zeros written
    
    //loop through rows
    for (int row = 0; row < n; row++){
        //allocate ja,ma for that row
        spmt.ja[row] = new int[spmt.nzcount[row]];
        spmt.ma[row] = new double[spmt.nzcount[row]];
        
        //copy from vectors
        std::copy(vec_ma.begin() + nzs,vec_ma.begin() + nzs + spmt.nzcount[row],spmt.ma[row]);
        std::copy(vec_ja.begin() + nzs,vec_ja.begin() + nzs +  spmt.nzcount[row],spmt.ja[row]);
        nzs += spmt.nzcount[row];
    }
}


//Convert a Mat into a sparse matrix;
void convert_to_CSR(const Mat &mat, SparMat &spmt){
    vector<int> vec_ja, vec_nzcount;
    vector<double> vec_ma;
    int n = mat.row_size;
    setupCS(&spmt,n,1);//allocate memory
    
    //build the appropriate vectors from Mat;
    int i = 0;
    int row_count = 0;
    float temp;
    
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
	vector<double> vec_ma;
	int n = gmap.size();
	setupCS(&spmt, n, 1);//allocate memory

	//build the appropriate vectors from Mat;
	for (auto node : gmap) {

		set<int> tempset = node.second; //adiacency list for the node.

		
		vec_ja.insert(vec_ja.end(), tempset.begin(), tempset.end()); //columns of nonzero elements = names of node children

		vector<double> tempvec(tempset.size(),1.); 
		vec_ma.insert(vec_ma.end(), tempvec.begin(), tempvec.end()); //entries = 1s. Graphmap are unweighted for now.
		
		vec_nzcount.push_back(tempset.size());//nonzero per row = number of node children
	}
	
	//use vectors to fill CSR
	fill_CSR(spmt, n, vec_nzcount, vec_ja, vec_ma);
}

void convert_to_MKL(SparMat &spmt, sparse_matrix_t &A){
    
    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
    
    int n = spmt.n;
    MKL_INT cols = n;
    MKL_INT rows = n;
    MKL_INT* rows_start= new MKL_INT[n];
    MKL_INT* rows_end = new MKL_INT[n];
    
    rows_start[0] = 0;
    rows_end[0] = spmt.nzcount[0];
    for (int i = 0; i<n;i++){
        rows_start[i] = spmt.nzcount[i-1] + rows_start[i - 1];
        rows_end[i] = spmt.nzcount[i] + rows_end[i - 1];
    }

    
    int nzs = rows_end[n-1];
    float* values = new float[nzs];
    MKL_INT* col_indx = new MKL_INT[nzs];
    
    int j = 0;
    for (int row = 0; row < n; row++){
        for (int i = 0; i <spmt.nzcount[row]; i++){
            col_indx[j] = spmt.ja[row][i];
            j++;
        }
    }
    
    j = 0;
    for (int row = 0; row < n; row++){
        for (int i = 0; i <spmt.nzcount[row]; i++){
            values[j] = spmt.ma[row][i];
            j++;
        }
    }
    
    
    mkl_sparse_s_create_csr (&A, indexing, rows, cols, rows_start,  rows_end, col_indx, values);
}

void permute(SparMat &spmt, int* perm){
    //permutes a matrix in CSR form
    int n = spmt.n;
    permute(spmt.nzcount,perm,n);
    permute(spmt.ja,perm,n);
    permute(spmt.ma,perm,n);
}

void matprint(const SparMat &spmt){
    cout << "PRINTING THE CSR MATRIX" << endl;
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
        }
        cout << mat.vec[i] << " ";
    }
    cout <<endl<< "finished printing the matrix"<<endl;

}

void matprint(vbsptr vbmat){
    int N = vbmat->n, *bsz = vbmat->bsz;
    int nBsj,sz,dim,col;
    cout << "PRINTING a VB MATRIX" << endl;
    cout<<"N ="<< N <<endl;
    for(int i = 0; i < N; i++ ) {
        cout<<"bsz["<<i+1<<"] = "<<bsz[i+1]<<endl;
        dim = bsz[i+1] - bsz[i];
        cout<<"nzcount["<<i<<"] = "<<vbmat->nzcount[i]<<endl;
        for(int j = 0; j<vbmat->nzcount[i]; j++){
            col = vbmat->ja[i][j];
            cout<<"ja["<<i<<"]["<<j<<"] = "<<col<<endl;
            nBsj = bsz[col];
            sz = bsz[col+1] - bsz[col];
            cout<<"BLOCK size="<<dim<<"x"<<sz<<endl;
            for(int k = 0; k < sz*dim; k++){
                if (k%dim != 0) cout<<" ";
                else cout<<endl;
                cout<<vbmat->ba[i][j][k];

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
	vector<double> vec_ma;
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
    vector<double> temp_mat(num_row*num_row, 0.0);
    
    // fill the matrix with data
    for (int l = 0; l < num_lines; l++)
    {
        double data;
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
    //loop vertically through blocks
    for(int i = 0; i < N; i++ ) {
        Hsz = bsz[i+1] - bsz[i];
        Bnum += vbmat->nzcount[i];
        //loop horizontaly through blocks
        for(int j = 0; j<vbmat->nzcount[i]; j++){
            col = vbmat->ja[i][j];
            Lsz = bsz[col+1] - bsz[col];
            BLsize.push_back(Lsz);
            BHsize.push_back(Hsz);
			tempSparse = 0;
            //loop through elements in the block
            for(int k = 0; k < Lsz*Hsz; k++){
                if (vbmat->ba[i][j][k] != 0.){
                    tempSparse += 1;
                }
            }
            Bsparsity.push_back(tempSparse);
        }
    }
}

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
int make_sparse_blocks(SparMat &spmt, VBSparMat &vbmat,double eps){
    int nBlock;
    int *nB = NULL, *perm = NULL;
    double *t_hash = NULL, *t_angle = NULL;
    
    if (init_blocks(&spmt, &nBlock, &nB, &perm, eps, t_hash, t_angle) != 0) {
        cout << "ERROR: COULD NOT CREATE PERMUTATION. Try with another epsilon" << endl;
        return -1;
    }
    
    permute(spmt, perm);

    int ierr = csrvbsrC(1, nBlock, nB, &spmt, &vbmat);
    if (ierr != 0){
        cout << "error occurred while creating block matrix" << endl;
        return -1;
    }
}

template <class myType>
void arrprint(myType *arr, int len){
    for(int i = 0; i < len; i++){
        cout<<arr[i]<<" ";
    }
    cout<<endl;
}
