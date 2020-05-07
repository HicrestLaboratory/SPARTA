#include "sparse_utilities.h"
#include "comp_mats.h"

typedef std::map<int, set<int> > Graphmap;


//GRAPH UTILITIES
//TODO RMAT reader/generator (maybe in Python)

void read_snap_format(Graphmap & gmap, string filename)
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

			if (gmap.count(current_node) == 0) { //new source node encountered
				gmap[current_node] = emptyset;
			}

			//get target node (column)
			getline(infile, temp);
			child = stoi(temp);
			gmap[current_node].insert(child);
		}
	}

//check if a graphmap is well-numbered, complete and (optional) symmetric
int isProper(const Graphmap & gmap, bool mirroring) {
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
void MakeUndirected(Graphmap & gmap) {

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
void MakeProper(Graphmap & gmap) {

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
void write_snap_format(Graphmap & gmap, string filename) {
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

void convert_to_mat(const Graphmap& gmap, DataT* mat, int mat_fmt) {

	int n = gmap.size();
	for (int i = 0; i < n; i++)
	{

	}



	vector<int> vec_ja, vec_nzcount;
	vector<DataT> vec_ma;
	int n = gmap.size();

	//build the appropriate vectors from Mat;
	for (auto node : gmap) {

		set<int> tempset = node.second; //adiacency list for the node.


		vec_ja.insert(vec_ja.end(), tempset.begin(), tempset.end()); //columns of nonzero elements = names of node children

		vector<DataT> tempvec(tempset.size(), 1.);
		vec_ma.insert(vec_ma.end(), tempvec.begin(), tempvec.end()); //entries = 1s. Graphmap are unweighted for now.

		vec_nzcount.push_back(tempset.size());//nonzero per row = number of node children
	}

}


	/*-------------------------------------------------------------------------------------------------------------
	END OF GRAPH UTILITIES
	*/