
#pragma once 

#include <sdsl/bit_vectors.hpp>
#include <cmph.h> 

#include<string>
using namespace std;
using namespace sdsl;


uint64_t CATEGORY_RUN=(uint64_t) 3;
uint64_t CATEGORY_COLCLASS=(uint64_t) 0;
uint64_t CATEGORY_COLVEC=(uint64_t) 2;
uint64_t CATEGORY_COLVEC_ONE = (uint64_t) 4; //100
uint64_t CATEGORY_COLVEC_TWO = (uint64_t) 5; //101


class OutputFile{
	public:
		string filename;
		std::ofstream fs;
	OutputFile(){

	}
	OutputFile(string filename){
		this->filename = filename;
		fs.open (filename.c_str(),  std::fstream::out );
	}
	void init(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::out);
	}
	void write(string towrite){
		fs << towrite; // <<endl;
	}
	void close(){
		fs.close();
	}
	~OutputFile(){
		fs.close();
	}
};
class DebugFile : public OutputFile	//derived class
{
	public:
		DebugFile(string filename){
			//if(!DEBUG_MODE){
					this->filename = filename;
					fs.open (filename.c_str(),  std::fstream::out );
			//
		}
		DebugFile(){
		}
};
class LogFile : public OutputFile	//derived class
{
	public:
		void log(string param_name, string param_value, string delim=":")
		{
			fs << param_name << delim << param_value << endl;
		}


};
class InputFile{
	public:
	string filename;
	std::fstream fs;
	InputFile(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::in);
	}
	InputFile(){
	}
	void init(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::in);
	}
	void rewind(){
		this->fs.close();
		this->fs.open(this->filename, fstream::in);
	}

	void close(){
		fs.close();
	}
	~InputFile(){
		fs.close();
	}
};
class Hashtable {
	public:
    std::unordered_map<uint32_t, uint32_t> htmap; // m_to_l
	uint32_t curr_id = 0;
	Hashtable(){
		curr_id = 0;
	}

	void copyFrom(Hashtable & h){
		if(htmap.size()!=0){
			htmap.clear();
		}
		for(const auto& entry:  h.htmap)
		{
			this->htmap[entry.first] = entry.second;
		}
		this->curr_id = h.curr_id;
	}

    uint64_t put_and_getid(uint32_t key) {
		if(htmap.count(key) > 0){ // present
			return htmap[key];
		}  else {	// absent
			htmap[key] = curr_id;
			curr_id+=1;
			return curr_id-1; 
		}
    }

	bool exists(uint32_t key){
		return htmap.count(key) > 0;
	}

	void clear(){
		if (htmap.size()!=0)
		{
			htmap.clear();
		}
		curr_id = 0;
	}

	vector<uint32_t> get_array(){
		vector<uint32_t> array(curr_id, 0);
		for (const auto& x : htmap){
			array[x.second] =  x.first  ;
		}
		return array;
	}
};

namespace CMPH{
	cmph_t *hash_cmph = NULL;
	void create_table(string filename, int nelem ){
		FILE * keys_fd = fopen(filename.c_str(), "r");
		
		if (keys_fd == NULL) 
		{
		fprintf(stderr, "File not found\n");
		exit(1);
		}	
		// Source of keys
		cmph_io_adapter_t *source = cmph_io_nlfile_adapter(keys_fd);
	
		cmph_config_t *config = cmph_config_new(source);
		cmph_config_set_algo(config, CMPH_CHM);
		hash_cmph = cmph_new(config);
		cmph_config_destroy(config);
		
		cmph_io_nlfile_adapter_destroy(source);   
		fclose(keys_fd);
	}

	unsigned int lookup(string str){	
		const char *key = str.c_str(); 
		//Find key
		unsigned int id = cmph_search(hash_cmph, key, (cmph_uint32)strlen(key));
		// fprintf(stderr, "Id:%u\n", id);
		//Destroy hash
		//cmph_destroy(hash);
		return id;
	}

	void mphf_destroy(){
		cmph_destroy(hash_cmph);
	}
}
using namespace CMPH;

typedef std::vector<bool> HuffCode;
typedef std::map<u_int32_t, HuffCode> HuffCodeMap;

namespace Huffman{
	/// @brief source rosetta code
	class INode
	{
		public:
			const int f;
			virtual ~INode() {}
		protected:
			INode(int f) : f(f) {}
	};
	class InternalNode : public INode
	{
	public:
		INode *const left;
		INode *const right;

		InternalNode(INode* c0, INode* c1) : INode(c0->f + c1->f), left(c0), right(c1) {}
		~InternalNode()
		{
			delete left;
			delete right;
		}
	};
	class LeafNode : public INode
	{
	public:
		u_int32_t c;

		LeafNode(int f, u_int32_t c) : INode(f), c(c) {}
	};

	struct NodeCmp
	{
		bool operator()(const INode* lhs, const INode* rhs) const { return lhs->f > rhs->f; }
	};

	INode* BuildTree(u_int32_t* frequencies, u_int32_t UniqueSymbols)
	{
		std::priority_queue<INode*, std::vector<INode*>, NodeCmp> trees;
	
		for (u_int32_t i = 0; i < UniqueSymbols; ++i)
		{
			if(frequencies[i] != 0)
				trees.push(new LeafNode(frequencies[i], (u_int32_t)i));
		}
		while (trees.size() > 1)
		{
			INode* childR = trees.top();
			trees.pop();

			INode* childL = trees.top();
			trees.pop();

			INode* parent = new InternalNode(childR, childL);
			trees.push(parent);
		}
		return trees.top();
	}

	void GenerateCodes(const INode* node, const HuffCode& prefix, HuffCodeMap& outCodes)
	{
		if (const LeafNode* lf = dynamic_cast<const LeafNode*>(node))
		{
			outCodes[lf->c] = prefix;
		}
		else if (const InternalNode* in = dynamic_cast<const InternalNode*>(node))
		{
			HuffCode leftPrefix = prefix;
			leftPrefix.push_back(false);
			GenerateCodes(in->left, leftPrefix, outCodes);

			HuffCode rightPrefix = prefix;
			rightPrefix.push_back(true);
			GenerateCodes(in->right, rightPrefix, outCodes);
		}
	}
}

using namespace Huffman;
//using namespace Helper;


struct simplitig {
	int nKmer; // number of k-mers in a simplitig
	int l; //local id: number of unique classes
	int optimal_bigD;  // optimal big D, delta value 0,1,2
	int optimal_useLocal;  // optimally use localTable if 1, 0 otherwise
	int optimal_space; //optimal space usage when using optimal D and useLocal
	bool mono = false;
	int length;
};

//
namespace TimeMeasure
{
	double time_start(){
        struct timeval timet;
        double t_begin;
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
        return t_begin;
    }
	double time_end(double t_begin, string msg){
        struct timeval timet;
        double t_end;
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
		if(msg!=""){
			cout<<msg<<" time = ";
			printf("%.2fs\n",t_end - t_begin);
		}
		return t_end;
	}
} 


namespace Helper{
	/**
	 * Compute mean for a collection of data.
	 * @param[in] v the data where T: float/int/unsigned int
	 * @return mean of `values`, or 0.0 if `values` is empty.
 	*/
	template <typename T> float get_average(vector<T> v){
		if(v.size()==0){
			return 0;
		}
		T summ = 0;
		for (T e:  v){
			summ+=e;
		}
		return summ/1.0/v.size();
	}

	void write_number_at_loc_advanced_by_block_sz(vector<uint64_t> & positions, uint64_t num, uint64_t loc_advanced_by_block_sz, uint64_t block_sz){ //requires loc_advanced_by_block_sz += block_size; 
		int64_t j=0;
		uint64_t begin = loc_advanced_by_block_sz;
		stack<uint64_t> qpositions;
		if(num!=0){
			if(block_sz==0)
				cout<<"must be non zero block size"<<endl;
		}
		while(num!=0)
		{
			if(num%2 == 1){
				//positions.push_back(loc_advanced_by_block_sz-1-j); //b[loc_advanced_by_block_sz-1+j] = num%2;
				qpositions.push(loc_advanced_by_block_sz-1-j);
			}
			j++;
			num /= 2;
		}
		while(!qpositions.empty()){
			positions.push_back(qpositions.top());
			qpositions.pop();
		}

		// if(DEBUG_MODE) debug1.fs<<-j<<" "<<block_sz<<endl;
		if (j > block_sz){
			cout<<"error in block"<<endl;
		}

	}

	void write_one(vector<uint64_t> & positions, uint64_t& b_it ){
		positions.push_back(b_it);
		b_it+=1;
	}

	void write_zero(vector<uint64_t> & positions, uint64_t& b_it ){
		b_it+=1;
	}
	
	void write_number_at_loc(vector<uint64_t> & positions, uint64_t num, uint64_t block_size, uint64_t& b_it ){
		write_number_at_loc_advanced_by_block_sz(positions, num, b_it+block_size, block_size);
		b_it += block_size; //successfully written and place on next bit; if size is 2, (0,1) written, now val is 2.
	
	}

	void write_colorbv_at_loc(vector<uint64_t> & positions, uint64_t& b_it, uint64_t* arr, int C ){
		 int nBlocks = ceil(C/64.0);
		 for(int i = 0; i < nBlocks; i++ ){
			write_number_at_loc(positions, arr[i], min(64, C), b_it );
			if(i==nBlocks-1){
				write_number_at_loc(positions, arr[i], C-64*i, b_it );
			}
		 }
	}

	void write_unary_one_at_loc(vector<uint64_t> & positions, uint64_t unary_num, uint64_t& b_it ){
		for(uint64_t i = 0; i<unary_num; i++ ){
			positions.push_back(b_it);
			b_it+=1;
		}
	}
	void write_unary_zero_at_loc(vector<uint64_t> & positions, uint64_t unary_num, uint64_t& b_it ){
		b_it+=unary_num;
	}

	void write_category(vector<uint64_t> & positions, uint64_t & b_it, uint64_t category, int bigD, int hd){ //0 run, case_lm, case_dlc
		if(category == CATEGORY_RUN){
			if(bigD==0){ //if(false){ //
				write_one(positions, b_it);
			}else{
				write_number_at_loc(positions, CATEGORY_RUN, (uint64_t)2, b_it); //11
			}
		}else if(category == CATEGORY_COLCLASS){
			if(bigD==0){
				write_zero(positions, b_it);
			}else{
				write_number_at_loc(positions, CATEGORY_COLCLASS, (uint64_t)1, b_it); //0
			}
		}else if(category == CATEGORY_COLVEC){ //assert bigD == 1
			if(bigD == 1 && hd == 1){
				write_number_at_loc(positions, CATEGORY_COLVEC, (uint64_t)2, b_it); //10
			}else if(bigD == 2 && hd==2){
				write_number_at_loc(positions, CATEGORY_COLVEC_TWO, (uint64_t)3, b_it);
			}else if(bigD == 2 && hd==1){
				write_number_at_loc(positions, CATEGORY_COLVEC_ONE, (uint64_t)3, b_it);
			}
		}
	}

	void write_binary_string_at_loc(vector<uint64_t> & positions, string binarystring, uint64_t& b_it){
		for (size_t i = 0; i< binarystring.length(); i++) {
			if (binarystring[i]=='1'){
				positions.push_back(b_it+i);
			}
		}
		b_it += binarystring.length();
	}

	void write_binary_vector_at_loc(vector<uint64_t> & positions, vector<bool> binary_vector, uint64_t& b_it){
		for (size_t i = 0; i< binary_vector.size(); i++) {
			if (binary_vector[i]== 1){
				positions.push_back(b_it+i);
			}
		}
		b_it += binary_vector.size();
	}

	void store_as_sdsl(vector<uint64_t>& positions, uint64_t bv_size, string filename){
		
		//bit_vector bv = bit_vector(bv_size, 0);
		bit_vector bv(bv_size, 0);
		uint64_t lastp  =0;
		for (uint64_t p: positions){
			bv[p] = 1;
		}
		if(filename=="rrr_main"){
			//if(DEBUG_MODE) debug2.fs<<bv;
		}
		rrr_vector<256> rrr_bv(bv);
		//cout << "rrr_MB_bv_mapping="<<size_in_bytes(rrr_bv_mapping)/1024.0/1024.0 << endl;
		store_to_file(rrr_bv, filename);	//"rrr_bv_mapping.sdsl"

		//return bv;
	}

	void store_as_binarystring(vector<uint64_t>& positions, uint64_t bv_size, string filename){
		OutputFile binarystring_file(filename);
		//sort(positions.begin(), positions.end());
		uint64_t bvi = 0;
		for (uint64_t k = 0; k<bv_size; k++){
			if(bvi < positions.size()){
				if(positions[bvi]==k){
					binarystring_file.write("1");
					bvi++;
				}else{
					binarystring_file.write("0");
				}
			}else{
				binarystring_file.write("0");
			}	
		}
	}

	
}
