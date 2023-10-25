#define VERSION_NAME "October 2023, Restructuring code 2"
#include <cmph.h> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <vector>
#include <stack>
#include <deque>
#include <fstream>
#include <algorithm>
#include <sdsl/bit_vectors.hpp>
#include<unordered_map>
#include<sstream>
#include<string>
#include <queue>
#include <map>
#include <climits> // for u_int32_t_BIT
#include <iterator>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "helper.hpp"
#include <unordered_map>


using namespace std;
using namespace sdsl;
using namespace TimeMeasure;
using namespace Helper;

const bool USE_TEST_METHOD = false;
bool NEW_DEBUG_MODE = false;
int DEBUG_MODE = 2;
const bool USE_MONO_BIT = true;

class COLESS{
public:
	float tot_time_hd=0;
	float tot_time_read=0;
	float tot_time_stoull_read=0;
	uint64_t sum_of_opt_space = 0;
	InputFile dedup_bitmatrix_file, spss_boundary_file, dup_bitmatrix_file;
	DebugFile debug1;
	DebugFile debug2;
	DebugFile all_ls;
	uint64_t num_kmers;
	int num_blocks_req;
	int num_simplitig;
	int M;
	int C;
	int max_run = 16;
	int lmaxrun;
	int max_run_choice = 1;
	vector<uint64_t> positions;
	HuffCodeMap huff_code_map;

	int lm = 0;
	int lc = 0;

	vector <struct simplitig> simplitigs; //l, optimal_bigD, optimal_useLocal, optimal_space 
	vector<char> spss_boundary; 

	//run param
	int d_class_diff = 1; //0,1,2
	bool USE_LOCAL_TABLE = true;
    bool USE_HUFFMAN = true;
	bool ALWAYS_LOCAL_OR_GLOBAL = false;

	COLESS(uint64_t num_kmers, int M, int C, string dedup_bitmatrix_fname, string dup_bitmatrix_fname, string spss_boundary_fname, int max_run){
		dedup_bitmatrix_file.init(dedup_bitmatrix_fname);
		spss_boundary_file.init(spss_boundary_fname);
		dup_bitmatrix_file.init(dup_bitmatrix_fname);
		this->max_run = max_run;
		this->num_kmers = num_kmers;
		this->M = M;
		this->C = C;
		this->lm = ceil(log2(M));
		this->lc = ceil(log2(C));
		this->num_blocks_req = ceil(C/64.0);
		num_simplitig = 0;
		if(DEBUG_MODE) debug1.init("debug1");
		if(DEBUG_MODE) debug2.init("debug2");
		if(DEBUG_MODE) all_ls.init("all_ls");
	}

	~COLESS(){
		mphf_destroy();
		// delete global_table;
	}

	inline void string_to_uint_colorbv(string line_of_bits, int nColor, uint64_t& lo, uint64_t& hi){ //lo means 0 [first 64 characters cut -c1-64], hi means 1
		//can also get nColor from line_of_bits
		double ts = time_start();
		lo = std::stoull(line_of_bits.substr(0,std::min(64,int(nColor))), nullptr, 2) ; 
		hi = 0;
		if(nColor > 64){
			string ss = line_of_bits.substr(64,(nColor-64));
			hi  = std::stoull(ss, nullptr, 2);
		}
		tot_time_stoull_read+=time_end(ts, "");
	}

	/// @brief WARNING: must be FREED
	void string_to_uint_colorbv(string line_of_bits, int nColor, uint64_t* bvarr){ //lo means 0 [first 64 characters cut -c1-64], hi means 1
		double ss = time_start();
		int num_blocks_req = ceil(nColor/64.0);
		//bvarr = new uint64_t[num_blocks_req];
		// little endian (lo=0,hi=num_blocks_req-1)
		bvarr[0] = std::stoull(line_of_bits.substr(0,std::min(64,int(nColor))), nullptr, 2) ;  //lo
		for(int i = 1; i<num_blocks_req; i++){
			bvarr[i] = 0;
			//range:  >= 64*i and < 64*i + 64
			if(i == num_blocks_req - 1){ 
				int rem = (nColor%64==0)?64:nColor%64; //if fractional vs whole
				string ss = line_of_bits.substr(64*i,rem);
				bvarr[i]   = std::stoull(ss, nullptr, 2);
			}else{
				string ss = line_of_bits.substr(64*i,64);
				bvarr[i]   = std::stoull(ss, nullptr, 2);
			}
		}
		tot_time_stoull_read+=time_end(ss, "");
	}

	int hammingDistance (uint64_t x, uint64_t y) {
		uint64_t res = x ^ y;
		return __builtin_popcountll (res) ;
	}

	/// @brief array version
	int hammingDistance (uint64_t*& x_bv_arr,uint64_t*& y_bv_arr, int num_blocks) {
		double ss = time_start();
		int sum  = 0;
		if(x_bv_arr==NULL){
			for(int i=0; i<num_blocks; i++){
				//#pragma vector always
				sum+= __builtin_popcountll (y_bv_arr[i]);
			}
			return sum;
		}else if(y_bv_arr==NULL){
			for(int i=0; i<num_blocks; i++){
				//#pragma vector always
				sum+= __builtin_popcountll (x_bv_arr[i]);
			}
			return sum;
		}else{
			for(int i=0; i<num_blocks; i++){
				//#pragma vector always
				sum+= hammingDistance(x_bv_arr[i], y_bv_arr[i]);
			}
		}
		tot_time_hd+=time_end(ss,"");
		return sum;
	}
	
	void copy_bv_arr(uint64_t*& from, uint64_t*& to, int num_blocks){
		if(to==NULL){
			to = new uint64_t[num_blocks];
		}
		for(int bid=0; bid<num_blocks; bid++){
			to[bid] = from[bid];
		}
	}

	inline uint64_t convert_block_ibit(int i_bit, int blockId, int C){
		int BYTESIZE=64;
		int num_blocks_req = ceil(C/64.0);
		if(blockId == num_blocks_req-1){ //last_block
			int rem = C%BYTESIZE;
			return blockId*BYTESIZE + (rem-1 - i_bit); //in range [0,C-1]
		}else{
			return  blockId*BYTESIZE + (BYTESIZE-1 - i_bit); //in range [0,C-1]
		}
	}
	/// @brief Do a sorted delta. gives a global table with slightly smaller size than just storing all color vectors
	void store_NONMST_global(){
		ifstream uniq_ms(dedup_bitmatrix_file.filename);
		string hd_line;
		
		vector<uint64_t> positions;
		uint64_t b_it = 0;

		vector<uint64_t> positions_hd;
		uint64_t b_it_hd = 0;

		uint64_t zero_64bit = 0;
		uint64_t prev_lo = 0;
		uint64_t prev_hi = 0;
		uint64_t* prev_bv_arr = NULL;
		uint64_t* bv_arr = new uint64_t[num_blocks_req];
		uint64_t lo, hi;

		//global_table = new string[M];
		int hdsum = 0;
		for(int i = 0 ; i< M; i++){
			string uniq_ms_line;
			getline(uniq_ms, uniq_ms_line);
			////
			////
			//unsigned int idx = lookup(uniq_ms_line);		// returns an if in range (0 to M-1) 
			//assert(idx < M);
			//global_table[idx] = uniq_ms_line;
			//assert(x==idx);
			////
			////
			//string_to_uint_colorbv(uniq_ms_line, C, lo, hi);
			//readBitVectorFromString(uniq_ms_line, bv);

			// read in C length value in bv array
			string_to_uint_colorbv(uniq_ms_line, C, bv_arr);
			if(DEBUG_MODE==3)cout<<uniq_ms_line<<"line is"<<endl;
			for(int blockId = 0; blockId<ceil(C/64.0); blockId++){
				if(DEBUG_MODE==3) cout<<bv_arr[blockId]<<endl;
			}
			if(true){ //i!=0
				int hd_prev = hammingDistance(bv_arr, prev_bv_arr, ceil(C/64.0));
				//int hd_prev = hammingDistance(hi, prev_hi) + hammingDistance(lo, prev_lo);
				hdsum += hd_prev;
				int lc = ceil(log2(C));

				// DOUBLE-CHECK
				for(int blockId = 0; blockId<ceil(C/64.0); blockId++){
					int howmanybits = 64;
					if(blockId==num_blocks_req-1){
						howmanybits = ((C%64)==0)?64:(C%64);
					}
					for (int i_bit = 0; i_bit < howmanybits; i_bit += 1)
					{
						if(prev_bv_arr == NULL){
							if ((bv_arr[blockId] >> i_bit) & 1)
							{
								write_number_at_loc(positions_hd, convert_block_ibit(i_bit,blockId,C), lc, b_it_hd); // i_bit is the different bit loc
								if(DEBUG_MODE==3) cout<<i_bit<<"blockId"<<blockId<<" "<<b_it_hd<<endl;
							}
						}else{
							if (((prev_bv_arr[blockId] >> i_bit) & 1) != ((bv_arr[blockId] >> i_bit) & 1))
							{
								write_number_at_loc(positions_hd, convert_block_ibit(i_bit,blockId,C), lc, b_it_hd); // i_bit is the different bit loc
								if(DEBUG_MODE==3) cout<<i_bit<<"blockId"<<blockId<<" "<<b_it_hd<<endl;
							}
						}	
						
					}
				}
				
				/*
				for (int i_bit = 0; i_bit < 64 && i_bit < C; i_bit += 1)
				{
					if (((prev_lo >> i_bit) & 1) != ((lo >> i_bit) & 1))
					{       
						write_number_at_loc(positions_hd, i_bit, lc, b_it_hd); // i_bit is the different bit loc
					}
				}
				for (int i_bit = 64; i_bit < C; i_bit += 1)
				{
					int actual_i_bit = i_bit - 64;
					if (((prev_hi >> actual_i_bit) & 1) != ((hi >> actual_i_bit) & 1))
					{
						write_number_at_loc(positions_hd, i_bit, lc, b_it_hd); // i_bit is the different bit loc
					}
				}
				*/

				// for(int ii = 0; ii< (hd_prev*lc)-1; ii++){
				//     write_zero(positions, b_it); // i_bit is the different bit loc
				// }
				// write_one(positions, b_it); // i_bit is the different bit loc
				
				if(hd_prev!=0){
					b_it+=hd_prev*lc-1; //skip these many locations: they are all 0
					if(DEBUG_MODE==3) cout<<"biter "<<b_it<<" "<<hd_prev<<endl;
					write_one(positions, b_it);
				}
				//if(i!=0) g.addEdge(i, i-1, hd_prev);
			}
			//prev_lo=lo;
			//prev_hi=hi;
			copy_bv_arr(bv_arr, prev_bv_arr, num_blocks_req);
		}
		uniq_ms.close();

		if(prev_bv_arr!=NULL) delete[] prev_bv_arr;
		if(bv_arr!=NULL) delete[] bv_arr;

		store_as_sdsl(positions_hd, b_it_hd, "rrr_map_hd" );	
		//write_binary_bv_from_pos_vector( positions_hd, b_it_hd, "rrr_map_hd" );
		store_as_sdsl(positions, b_it, "rrr_map_hd_boundary" );	
		
	}

	void store_global_color_class_table(){
		vector<uint64_t> positions;  //wasteful
		uint64_t b_it = 0;

		//LogFile log_num_color_in_class;
		//log_num_color_in_class.init("log_num_color_in_class"); 
		dedup_bitmatrix_file.rewind();
		//global_table = new string[M];bm::
		for(int x=0; x<M; x++){
			string bv_line;
			getline(dedup_bitmatrix_file.fs, bv_line);
			unsigned int idx = lookup(bv_line);		// returns an if in range (0 to M-1) 
			assert(idx < M);
			//global_table[idx] = bv_line;
			assert(x==idx);

			uint64_t* colbv = new uint64_t[num_blocks_req];
			string_to_uint_colorbv(bv_line, C, colbv);
			//bm::bvector<> colbv;
			//readBitVectorFromString(bv_line, colbv);

			write_colorbv_at_loc(positions, b_it, colbv, C);

			// write_number_at_loc(positions, array_lo[idx], min(64, C), b_it ); //array_hi[x] higher uint64_t
			// write_number_at_loc(positions, array_hi[idx], C-64, b_it ); //array_lo[x] lower uint64_t
			delete[] colbv;
		}
		dedup_bitmatrix_file.fs.close();

		//store_as_binarystring(positions, b_it, "bb_map" );
		store_as_sdsl(positions, b_it, "rrr_map" );
		store_as_sdsl(positions, b_it, "rrr_bv_mapping.sdsl");
		cout << "expected_MB_bv_mapping="<<(C*M)/8.0/1024.0/1024.0 << endl;
		//cout << "rrr_MB_bv_mapping="<<sdsl::size_in_bytes(store_as_sdsl(positions, b_it, "rrr_bv_mapping.sdsl" ))/1024.0/1024.0 << endl;
	}

	bool startKmer(vector<char>& spss_boundary, uint64_t it_kmer){
		return spss_boundary[it_kmer] == '1';
	}

	bool endKmer(vector<char>& spss_boundary, uint64_t it_kmer){
		uint64_t num_kmers = spss_boundary.size();
		return spss_boundary[(it_kmer + 1) % num_kmers] =='1';
	}

	struct local_per_simplitig{
		vector<int> hds;
		vector<unsigned int> curr_kmer_cc_id;
		int space_needed;
	};

	void method1_pass1()
	{ 
		
		string combo_string = "";
		//DebugFile cases_smc("cases_smc");
		//DebugFile debuglll("lll");
		DebugFile skipper("skipper");
		DebugFile debug_combo;
		
		if(DEBUG_MODE) skipper.init("skipper");
		if(DEBUG_MODE) debug_combo.init("combo");

		double ts = time_start();
		create_table(dedup_bitmatrix_file.filename, M );
		time_end(ts, "CMPH constructed perfect hash for "+to_string(M)+" keys.");

		ts = time_start();
		bool skip_global_load=false;
		if(NEW_DEBUG_MODE==true){
			skip_global_load=true;
		}
		
		//OutputFile cmp_keys;
		if(skip_global_load==false){
			//cmp_keys.init("cmp_keys");  // get frequency count
		}

		
		vector<uint64_t> frequencies_colclass(M, 0);
		uint32_t prev_col_class, curr_col_class;

		//convert to color class from color vector
		for (uint64_t i=0; i < num_kmers; i+=1){  // read two files of length num_kmers 
			string spss_line, bv_line;
			getline (spss_boundary_file.fs,spss_line); 
			spss_boundary.push_back(spss_line[0]); //this kmer starts a simplitig
			
			if(skip_global_load==false) getline (dup_bitmatrix_file.fs,bv_line);
			//if(skip_global_load==false) cmp_keys.fs << lookup(bv_line) <<endl;
			curr_col_class = lookup(bv_line);
			if(spss_line[0]=='1'){	//start
				num_simplitig += 1;
				if(skip_global_load==false) frequencies_colclass[curr_col_class] +=1;
			}else{
				if(prev_col_class!=curr_col_class){
					if(skip_global_load==false) frequencies_colclass[curr_col_class] +=1;
				}
			}
			prev_col_class = curr_col_class;
		}


		vector<struct simplitig> temp_simplitigs(num_simplitig);
		simplitigs = temp_simplitigs;


		uint64_t prev_end = 0;
		int simplitig_it = 0;
		for (uint64_t i=0; i < num_kmers; i+=1){ 
			
			
			if(endKmer(spss_boundary, i)){
				if(prev_end==0) {
					simplitigs[simplitig_it++].length = i-prev_end+1;
				}else{
					simplitigs[simplitig_it++].length = i-prev_end;
				}
				
				//cout<<i<<" "<<prev_end<<endl;
				prev_end = i;
			}
		}
		time_end(ts, "CMPH lookup for "+to_string(num_kmers)+"keys.");
		//if(skip_global_load==false) cmp_keys.close();

		OutputFile outfile_freq("frequency_sorted");
		for (uint32_t i=0; i < M; i+=1){ 
			outfile_freq.fs<<frequencies_colclass[i]<<endl;
		}
		if(skip_global_load==false){
			double ts = time_start();
			//system("cat cmp_keys | sort -n | uniq -c | rev | cut -f 2 -d\" \" | rev > frequency_sorted");
			
			time_end(ts, "Sorting and getting freq for "+to_string(num_kmers)+" keys.");
		}
		
		outfile_freq.close();
		

		ts = time_start();
		InputFile infile_freq("frequency_sorted");
		string line;
		// Build frequency table
		u_int32_t *frequencies = new u_int32_t[M]; // M -> no. of unique symbols
		std::fill_n(frequencies, M, 0);
		u_int32_t i = 0;
		while(getline(infile_freq.fs, line)){
			stringstream ss(line);
			u_int32_t a ;
			ss >> a; 
			frequencies[i++]= a;
		}		
		time_end(ts, "Read freq for "+to_string(M)+" values.");
		infile_freq.close();

		ts = time_start();
		INode* root = BuildTree(frequencies, M);
        GenerateCodes(root, HuffCode(), huff_code_map); // huff_code_map is filled: uint32t colclassid-> vector bool
		delete frequencies;
		delete root;
		time_end(ts, "Build huffman tree on " +to_string(M)+" values.");

		ts = time_start();
		store_NONMST_global();
		//store_global_color_class_table();
		time_end(ts, "Written global table for "+to_string(M)+" values.");

		string bv_line;
		DebugFile optout;
		if(DEBUG_MODE) optout.init("optout");

		vector<uint32_t> optimal_ht;
		lmaxrun = ceil(log2(max_run));

		// per simplitig values
		Hashtable local_hash_table;
		uint64_t sum_length_huff_nonrun = 0;
		uint64_t sum_length_huff_uniq_nonrun = 0;
		uint64_t sum_dlc_space = 0;
		uint64_t sum_skip_space = 0;

		uint64_t skip = 0;
		int case_run = 0;
		int case_lm = 0;
		int case_dlc = 0;
		//
		vector<uint64_t> positions_local_table;
		vector<uint64_t> positions_mono;

		uint64_t b_it_local_table = 0;
		uint64_t b_it_mono = 0;

		// per kmer values
		uint64_t simplitig_start_id = 0;
		simplitig_it = 0;

		uint64_t prev_bv_lo = 0;
		uint64_t prev_bv_hi = 0;
		uint64_t* prev_bv = NULL;
		uint64_t* curr_bv = new uint64_t[num_blocks_req] ;
		
		struct local_per_simplitig local_simplitig;
		local_simplitig.hds.clear();
		local_simplitig.curr_kmer_cc_id.clear();


		for (int x = 0; x < num_simplitig; x++)
		{
			simplitigs[x].optimal_space = 99999999;
		}
		//	cout<<"U B C S "<<useLocal<<" "<<bigD<<" "<<curr_kmer_cc_id<<" "<<simplitig_it<<endl;		for (; big_d_local_combo < 6; big_d_local_combo++)
		int big_d_local_combo = 0;
		uint64_t it_kmer = 0;
		uint64_t local_it = 0;

		dup_bitmatrix_file.rewind();
		while (true)
		{
			//start with it-kmer 0
			// if (DEBUG_MODE)
			// 	all_ls.fs << "Start_bigd"
			// 			  << " " << big_d_local_combo <<it_kmer<<" "<<simplitig_it<<" " << endl;
			if( big_d_local_combo == 0 && startKmer(spss_boundary, it_kmer)){ //start of simp
				local_it = 0;
				local_simplitig.hds.clear();
				local_simplitig.curr_kmer_cc_id.clear();

				// just populate color class id for current simplitig
				while(true){
					getline (dup_bitmatrix_file.fs,bv_line);
			
					//string_to_uint_colorbv(bv_line, C, curr_bv_lo, curr_bv_hi);
					string_to_uint_colorbv(bv_line, C, curr_bv);
					// cout<<bv_line<<" bv line is"<<endl;
					// for(int blockId = 0; blockId<ceil(C/64.0); blockId++){
					// 	cout<<curr_bv[blockId]<<endl;
					// }

					if(not startKmer(spss_boundary, it_kmer+local_it)){ //non-start, could be end
						int hd = hammingDistance(prev_bv, curr_bv, ceil(C/64.0)) ;
						//cout<<hd<<endl;
						
						//int hd = hammingDistance(prev_bv_hi, curr_bv_hi) + hammingDistance(prev_bv_lo, curr_bv_lo);
						local_simplitig.hds.push_back(hd);
					}

					local_simplitig.curr_kmer_cc_id.push_back(lookup(bv_line)); // uint64_t num = bphf->lookup(curr_bv);

					// prev_bv_hi = curr_bv_hi;
					// prev_bv_lo = curr_bv_lo;
					//prev_bv = curr_bv;
					copy_bv_arr(curr_bv, prev_bv, num_blocks_req);
					
					if(endKmer(spss_boundary, it_kmer+local_it)){
						local_it = 0;
						break;
					}
						
					local_it++;
				}	
			}
			
			int useLocal = (big_d_local_combo / 3);
			int bigD = big_d_local_combo % 3;
			int hd = 0;
			unsigned int curr_kmer_cc_id = local_simplitig.curr_kmer_cc_id[local_it++];

			//cout << "U B C S " << useLocal << " " << bigD << " " << curr_kmer_cc_id << " " << simplitig_it << endl;
			if (not startKmer(spss_boundary, it_kmer))
			{ // non-start
				hd = local_simplitig.hds[local_it-2];
				if (hd == 0)
				{ // CAT=RUN
					skip += 1;
					case_run += 1;
				}
				else
				{ // CAT=NRUN

					if(skip!=0){
						//cout<<skip<<" "<<" "<<simplitigs[simplitig_it].length<<" "<<simplitig_it<<endl;
						skipper.fs<<skip<<endl;
						if(bigD==0){
							sum_skip_space += 1; 
						}else{
							sum_skip_space += 2; 
						}
						if(USE_TEST_METHOD){

							if(skip <= 4){
								max_run_choice = 0;
								max_run =  4;
							}else if(skip <= 16 ){
								max_run_choice = 1;
								max_run = 16;
							} else if (skip <= 128){
								max_run_choice = 2;
								max_run = 128;
							} else{
								max_run_choice = 3;
								max_run = 256;
							}
							lmaxrun = ceil(log2(max_run));
							sum_skip_space += 2;
						}
						sum_skip_space += floor(skip / max_run) + 1 + lmaxrun ;

					}
					skip = 0;

					if (hd <= bigD)
					{ // CAT=LC
						case_dlc += 1;
						if(bigD==2){
							sum_dlc_space += hd * lc + 3; // 101 = d>2 = 100 = d>1, d upto 4 per simp 2 bit
						}else{
							sum_dlc_space += hd * lc + 2; // 101 = d>2 = 100 = d>1, d upto 4 per simp 2 bit
						}
					}
					else
					{ // CAT=LM
						if (useLocal == 1)
						{
							local_hash_table.put_and_getid(curr_kmer_cc_id);
						}
						else
						{
							sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();
						}
						case_lm += 1;
					}
				}
			}
			else
			{ // start of simplitig, so CAT=LM
				simplitig_start_id = it_kmer;
				skip = 0;
				case_lm += 1;
				if (useLocal == 1)
				{
					local_hash_table.put_and_getid(curr_kmer_cc_id);
				}
				else
				{
					sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();
				}
			}

			if (endKmer(spss_boundary, it_kmer))
			{ // end k-mer of simplitig
				local_it = 0; 
				int l = -1;
				int ll = -1;

				if (useLocal == 1)
				{
					l = local_hash_table.curr_id;
					ll = ceil(log2(l) * 1.0);
					vector<uint32_t> local_ht_arr = local_hash_table.get_array();
					for (uint32_t i = 0; i < local_hash_table.curr_id; i++)
					{
						uint32_t uniq_col_class_id = local_ht_arr[i];
						sum_length_huff_uniq_nonrun += huff_code_map[uniq_col_class_id].size();
					}
				}

				if(skip!=0){
					skipper.fs<<skip<<endl;
					if(skip==simplitigs[simplitig_it].length-1){
							//cout<<skip<<" "<<simplitig_it<<" "<<simplitigs[simplitig_it].length<<endl;
							simplitigs[simplitig_it].mono = false;
					}
					if(bigD==0){
						sum_skip_space += 1; 
					}else{
						sum_skip_space += 2; 
					}
					if(USE_TEST_METHOD){
						if(skip <= 4){
							max_run_choice = 0;
							max_run =   4;
						}else if(skip <= 16 ){
							max_run_choice = 1;
							max_run = 16;
						} else if (skip <= 128){
							max_run_choice = 2;
							max_run = 128;
						} else{
							max_run_choice = 3;
							max_run = 256;
						}
						lmaxrun = ceil(log2(max_run));
						sum_skip_space += 2;
					}
					sum_skip_space += floor(skip / max_run) + 1 + lmaxrun ;
				}
				skip = 0;

				local_simplitig.space_needed = useLocal * (2 + 1 + lm + sum_length_huff_uniq_nonrun + (ll+1) * case_lm + sum_dlc_space  + sum_skip_space) + (1 - useLocal) * (2 + 1 + sum_length_huff_nonrun + sum_dlc_space + case_lm + sum_skip_space);
				
				if(DEBUG_MODE)
					optout.fs << "every: simp:"<<simplitig_it<<"bigD:"<< bigD<<" ul:"<<useLocal<<" space:"<<local_simplitig.space_needed<<" optbigD:"<< simplitigs[simplitig_it].optimal_bigD << " optLocal:" << simplitigs[simplitig_it].optimal_useLocal << " opspace:" << simplitigs[simplitig_it].optimal_space <<" sum_huff:"<<sum_length_huff_uniq_nonrun<<" sum_dlc: "<<sum_dlc_space<<"sum_skip_space: "<<sum_skip_space << endl;

				//if(per_simplitig_optimal_space[simplitig_it] == big_d_local_combo)//random
				
				if (local_simplitig.space_needed < simplitigs[simplitig_it].optimal_space)
				{
					if(useLocal==1){
						optimal_ht.clear();
						optimal_ht = local_hash_table.get_array();
						simplitigs[simplitig_it].l = local_hash_table.curr_id;
					}else{
						simplitigs[simplitig_it].l = 0;
					}
					simplitigs[simplitig_it].optimal_space = local_simplitig.space_needed;
					simplitigs[simplitig_it].optimal_bigD = bigD;
					simplitigs[simplitig_it].optimal_useLocal = useLocal;
				}
				
				case_run = case_lm = case_dlc = 0;
				sum_length_huff_nonrun = sum_length_huff_uniq_nonrun = sum_dlc_space = sum_skip_space = 0;

				if (big_d_local_combo < 5)
				{
					it_kmer = simplitig_start_id;
					local_it = 0;
					if(useLocal == 1)
						local_hash_table.clear();
					big_d_local_combo++;
					continue;
				}
				else
				{
					//it_kmer++;
					combo_string+=to_string(simplitig_it)+" "+to_string(simplitigs[simplitig_it].optimal_bigD)+" "+to_string(simplitigs[simplitig_it].optimal_useLocal)+"\n";
					if(simplitigs[simplitig_it].mono){
						////write_number_at_loc(positions_local_table,3, 2, b_it_local_table);
						write_one(positions_mono, b_it_mono);
					}else{
						write_number_at_loc(positions_local_table, simplitigs[simplitig_it].optimal_bigD, 2, b_it_local_table);
						write_zero(positions_mono, b_it_mono);
					}
					
					if(simplitigs[simplitig_it].mono){
						//write_one(positions_local_table, b_it_local_table);

					}else if (simplitigs[simplitig_it].optimal_useLocal == 1)
					{
						write_one(positions_local_table, b_it_local_table);
						write_number_at_loc(positions_local_table, optimal_ht.size(), lm, b_it_local_table);
						
						for (uint32_t ii = 0; ii < optimal_ht.size(); ii++)
						{
							uint32_t uniq_col_class_id = optimal_ht[ii];
							write_binary_vector_at_loc(positions_local_table, huff_code_map[uniq_col_class_id], b_it_local_table);
						}
						optimal_ht.clear();
					}
					else
					{
						write_zero(positions_local_table, b_it_local_table);
					}

					// if(DEBUG_MODE)
					// 	optout.fs << "curr: simp:"<<simplitig_it<<"bigD:"<< bigD<<" ul:"<<useLocal<< " optbigD:"<< per_simplitig_optimal_bigD[simplitig_it] << " optLocal:" << per_simplitig_optimal_useLocal[simplitig_it] << " opspace:" << per_simplitig_optimal_space[simplitig_it] << endl;

					// re-init for new simplitig
					if(useLocal==1)
						local_hash_table.clear();
					
					sum_of_opt_space+=simplitigs[simplitig_it].optimal_space;
					simplitig_it += 1;


					if (it_kmer != num_kmers)
						big_d_local_combo = 0;
				}
			}

			it_kmer++;			
			if (it_kmer == num_kmers)
				break;
		}
		if(prev_bv!=NULL) delete[] prev_bv;
		if(curr_bv!=NULL) delete[] curr_bv;

		cout << "b_it_local_table_size: " << b_it_local_table << endl;
		//store_as_binarystring(positions_local_table, b_it_local_table, "bb_local_table");
		store_as_sdsl(positions_local_table, b_it_local_table, "rrr_local_table");
		store_as_sdsl(positions_mono, b_it_mono, "rrr_mono");
		
		debug_combo.fs<<combo_string;
	}


	/// @brief writing phase
	void method1_pass2()
	{
		vector<uint64_t> positions;
		uint64_t b_it = 0;
		dup_bitmatrix_file.rewind();
		//DebugFile cases_smc("cases_smc");
		//DebugFile cases_skip("cases_skip");

		uint64_t curr_bv_hi = 0;
		uint64_t curr_bv_lo = 0;
		uint64_t prev_bv_hi = 0;
		uint64_t prev_bv_lo = 0;
		uint64_t* curr_bv = new uint64_t[num_blocks_req] ;
		uint64_t* prev_bv = NULL;

		//bm::bvector<> curr_bv;
		//bm::bvector<> prev_bv;

		uint64_t skip = 0;

		//InputFile cmp_keys("cmp_keys");
		uint64_t simplitig_it = 0;
		int l = simplitigs[0].l;
		int ll = ceil(log2(l));
		int lm_or_ll;
		
		Hashtable local_ht;
		lmaxrun = ceil(log2(max_run));
		for (uint64_t i = 0; i < num_kmers; i += 1)
		{
			int bigD = simplitigs[simplitig_it].optimal_bigD;
			int useLocal = simplitigs[simplitig_it].optimal_useLocal;
			l = simplitigs[simplitig_it].l;
			ll = ceil(log2(l));
			
			// if(DEBUG_MODE)
			// 	all_ls.fs<<bigD<<" "<<useLocal<<" "<<l<<" "<<ll<<endl;
			if (useLocal == 1)
			{
				lm_or_ll = ll;
			}
			else
			{
				lm_or_ll = lm;
			}
			
			// load the color vector of current k-mer from disk to "curr_bv_hi/lo"
			string bv_line;
			getline(dup_bitmatrix_file.fs, bv_line); // bv line = color vector C bits
			//string_to_uint_colorbv(bv_line, C, curr_bv_lo, curr_bv_hi);
			//readBitVectorFromString(bv_line, curr_bv);
			string_to_uint_colorbv(bv_line, C, curr_bv);
			// for(int blockId = 0; blockId<ceil(C/64.0); blockId++){
			// 			cout<<" ll:"<<curr_bv[blockId]<<endl;
			// 		}

			unsigned int curr_kmer_cc_id = lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
		
			if (not startKmer(spss_boundary, i))
			{ // non-start
				//int hd = hammingDistance(prev_bv_hi, curr_bv_hi) + hammingDistance(prev_bv_lo, curr_bv_lo);
				int hd = hammingDistance(prev_bv, curr_bv, ceil(C/64.0));
				//int hd = hd_bv(prev_bv, skip = curr_bv);

				if (hd == 0)
				{ // CATEGORY=RUN
					skip += 1;
					// case_run+=1;
					// if(DEBUG_MODE) cases_smc.fs << "r" << endl;
				}
				else
				{ // CATEGORY=NOT_RUN
					if (skip != 0)
					{ // not skipped, run break, write lm
						// paul method
						if(!simplitigs[simplitig_it].mono){
							if(USE_TEST_METHOD){
								// if(DEBUG_MODE) cases_skip.fs << skip << endl;
								if(skip <= 4){
									max_run_choice = 0;
									max_run =   4;
								}else if(skip <= 16 ){
									max_run_choice = 1;
									max_run = 16;
								} else if (skip <= 128){
									max_run_choice = 2;
									max_run = 128;
								} else{
									max_run_choice = 3;
									max_run = 256;
								}
								lmaxrun = ceil(log2(max_run));
								write_number_at_loc(positions, (uint64_t)max_run_choice, (uint64_t)2, b_it);

							}
							int q, rem;
							q = floor(skip / max_run);
							rem = skip % max_run;
							assert(skip == q * max_run + rem); // skip = q*max_run + rem
							write_category(positions, b_it, CATEGORY_RUN, bigD, 0);
							write_unary_one_at_loc(positions, (uint64_t)q, b_it);
							write_zero(positions, b_it);
							write_number_at_loc(positions, (uint64_t)rem, (uint64_t)lmaxrun, b_it);	
						}
					}
					skip = 0;

					if (hd <= bigD and !simplitigs[simplitig_it].mono)
					{ // CATEGORY=LC
						// if(hd*(lc + 1) < huff_code_map[curr_kmer_cc_id].size() && hd==1 ){ //CATEGORY=LC
						// if(hd*(lc + 1) < lm && hd==1){ //CATEGORY=LC
						// if(DEBUG_MODE) cases_smc.fs << "d" << endl;

						//category colvec = 100, 101 //cAT COLVEC = 10 HD = 1 
						//write_number_at_loc(positions, CATEGORY_COLVEC, 2, b_it);

						write_category(positions, b_it, CATEGORY_COLVEC, bigD, hd);
						
						/* block of bm::bvector*/
						/*
						unsigned i_bit = prev_bv.get_first(); //xor => prev_bv
						int count_bigD = 0;
						do
						{
							//cout << i_bit << endl;
							write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
							i_bit = prev_bv.get_next(i_bit);
							if (!i_bit ){
								break;
							}
						} while(1);
						*/
						for(int blockId = 0; blockId<ceil(C/64.0); blockId++){
							int howmanybits = 64;
							if(blockId==num_blocks_req-1){
								howmanybits = ((C%64)==0)?64:(C%64);
							}
							
							for (int i_bit = 0; i_bit < howmanybits; i_bit += 1)
							{
								if(prev_bv == NULL){
									if ((curr_bv[blockId] >> i_bit) & 1)
									{
										write_number_at_loc(positions, convert_block_ibit( i_bit,  blockId,  C), lc, b_it); // i_bit is the different bit loc
									}
								}else{
									if (((prev_bv[blockId] >> i_bit) & 1) != ((curr_bv[blockId] >> i_bit) & 1))
									{
										write_number_at_loc(positions, convert_block_ibit( i_bit,  blockId,  C), lc, b_it); // i_bit is the different bit loc
									}
								}	
								
							}
						}
				
						// for(int blockId = 0; blockId<ceil(C/64.0); blockId++){
						// 	int howmanybits = min(64, C);
						// 	// if(blockId == ceil(C/64.0)-1 ){
						// 	// 	howmanybits = C%64;
						// 	// }
						// 	for (int i_bit = 0; i_bit < howmanybits; i_bit += 1)
						// 	{
						// 		if (((prev_bv[i] >> i_bit) & 1) != ((curr_bv[i] >> i_bit) & 1))
						// 		{
						// 			write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
						// 		}
						// 	}
						// }
						

					// 	for (int i_bit = 0; i_bit < 64 && i_bit < C; i_bit += 1)
					// 	{
					// 		if (((prev_bv_lo >> i_bit) & 1) != ((curr_bv_lo >> i_bit) & 1))
					// 		{
					// 			write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
					// 		}
					// 	}
					// 	for (int i_bit = 64; i_bit < C; i_bit += 1)
					// 	{
					// 		int actual_i_bit = i_bit - 64;
					// 		if (((prev_bv_hi >> actual_i_bit) & 1) != ((curr_bv_hi >> actual_i_bit) & 1))
					// 		{
					// 			write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
					// 		}
					// 	}
					// 
					}else
					{ // CATEGORY=LM
						if(!simplitigs[simplitig_it].mono){
							write_category(positions, b_it, CATEGORY_COLCLASS, bigD, hd);
							if (useLocal==1)
							{
								// if(DEBUG_MODE) cases_smc.fs << "l" << endl;
								uint64_t localid = local_ht.put_and_getid(curr_kmer_cc_id);
								if (ll == 0 && localid == 1)
								{
									cout << "trouble" << endl;
								}
								write_number_at_loc(positions, localid, ll, b_it);
							}
							else
							{
								// if(DEBUG_MODE) cases_smc.fs << "m" << endl;
								if (USE_HUFFMAN)
								{
									write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);
								}
								else
								{
									write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
								}
							}
						}
					}
				}
			}
			else
			{ // start of simplitig, so CAT=LM
				l = simplitigs[simplitig_it].l;
				ll = ceil(log2(l));
				lm_or_ll = ll;

				// case_lm+=1;
				if(simplitigs[simplitig_it].mono){
					//cout<<simplitig_it<<" "<<l<<endl;
					write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);
					//continue;
				}else{
					//write_number_at_loc(positions, CATEGORY_COLCLASS, 1, b_it);
					write_category(positions, b_it, CATEGORY_COLCLASS, bigD, 0);
					if (useLocal == 1)
					{
						// if(DEBUG_MODE) cases_smc.fs << "l" << endl;
						uint64_t localid = local_ht.put_and_getid(curr_kmer_cc_id);
						if (ll == 0)
						{
							assert(localid == 0);
						}
						write_number_at_loc(positions, localid, ll, b_it);
					}
					else
					{
						// if(DEBUG_MODE) cases_smc.fs << "m" << endl;
						if (USE_HUFFMAN==true){
							write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);
						}else{
							write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
						}	
						// assert(curr_kmer_cc_id<M && curr_kmer_cc_id>0);
					}
				}
			}

			if (endKmer(spss_boundary, i))
			{ // end k-mer of simplitig
				local_ht.clear();
				simplitig_it += 1;
				if (useLocal){
					lm_or_ll = ll;
				}else{
					lm_or_ll = lm;
				}
				if (skip != 0  and !simplitigs[simplitig_it].mono)
				{ // not skipped, run break, write lm
					if(USE_TEST_METHOD){
						// if(DEBUG_MODE) cases_skip.fs << skip << endl;
						if(skip <= 4){
							max_run_choice = 0;
							max_run =   4;
						}else if(skip <= 16 ){
							max_run_choice = 1;
							max_run = 16;
						} else if (skip <= 128){
							max_run_choice = 2;
							max_run = 128;
						} else{
							max_run_choice = 3;
							max_run = 256;
						}
						lmaxrun = ceil(log2(max_run));
						write_number_at_loc(positions, (uint64_t)max_run_choice, (uint64_t)1, b_it);

					}
					int q, rem;
					q = floor(skip / max_run);
					rem = skip % max_run;
					assert(skip == q * max_run + rem); // skip = q*max_run + rem
					// paul method
					write_category(positions, b_it, CATEGORY_RUN, bigD, 0);
					write_unary_one_at_loc(positions, (uint64_t)q, b_it);
					write_zero(positions, b_it);
					write_number_at_loc(positions, (uint64_t)rem, (uint64_t)lmaxrun, b_it);
					// my method //100001
				}
				skip = 0;
			}
			prev_bv_hi = curr_bv_hi;
			prev_bv_lo = curr_bv_lo;
			copy_bv_arr(curr_bv, prev_bv, num_blocks_req);

		}
		if(prev_bv!=NULL) delete[] prev_bv;
		if(curr_bv!=NULL) delete[] curr_bv;

		// DebugFile positions_out("positions_out");
		// for (uint64_t tt : positions)
		// {
		// 	if(DEBUG_MODE) positions_out.fs << tt << endl;
		// }
		cout << "b_it_size: " << b_it << endl;
		store_as_sdsl(positions, b_it, "rrr_main");
		cout << "sum_of_opt_space: " << sum_of_opt_space << endl;
		
		//store_as_binarystring(positions, b_it, "bb_main");
	}
};

int main (int argc, char* argv[]){
	cout<<"Version: "<<VERSION_NAME<<endl;
	//srand(time(nullptr));
	//srand(0);
	vector<string> args(argv + 1, argv + argc);
    string dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname;
	//string tmp_dir;
    int M, C;
	int max_run = 16;
	
	uint64_t num_kmers=0;
	
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: tool -g debug -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run>" << endl;
            return 0;
        } else if (*i == "-i") {
            dedup_bitmatrix_fname = *++i;
        } else if (*i == "-d") {
            dup_bitmatrix_fname = *++i;
        }else if (*i == "-c") {
            C = std::stoi(*++i);
        }else if (*i == "-m") {
            M = std::stoi(*++i);
        }else if (*i == "-k") {
            num_kmers = std::stol(*++i);
        }else if (*i == "-s") {
            spss_boundary_fname = *++i;
		}else if (*i == "-x") {
            max_run = std::stoi(*++i);
		}else if (*i == "-g") {
            NEW_DEBUG_MODE = true;
		}
    }
	max_run = 16;
	COLESS coless(num_kmers, M, C, dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname, max_run);
	
	double ts = time_start();
	coless.method1_pass1();
	time_end(ts, "pass1.");
	cout<<"HD: "<<coless.tot_time_hd<<endl;
	cout<<"Read: "<<coless.tot_time_read<<endl;
	cout<<"Read stoull: "<<coless.tot_time_stoull_read<<endl;

	
	//exit(1);
	if(NEW_DEBUG_MODE==false){
		ts = time_start();
		coless.method1_pass2();
		time_end(ts, "pass2.");
	}

	return EXIT_SUCCESS;
}