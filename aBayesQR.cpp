#include <algorithm>
#include <random>
#include <iomanip> // std::setprecision
#include <vector> // for use of vector
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream> 
#include <numeric>
#include <time.h>


#define MAX_NUM_BRANCHES 10000
using namespace std;

int nChoosek(int n,int k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

int maxNforKcomb(int limit, int k)
{
    int n=k;
    while (nChoosek(n+1,k) < limit){
    n=n+1;
    }
    return n;

}

double delta0 = 0.1, delta1 = 0.001;

string binary(int number, stringstream& strs) 
{
	int remainder;
	if(number <= 1) 
	{
		strs << number;
	  	return strs.str() ;
	}
	remainder = number%2;
	binary(number >> 1, strs);    
	strs << remainder;
	
	return  strs.str();
}

// Code for parsing SAM file (line 59-339) adapted from PredictHaplo
// http://bmda.cs.unibas.ch/software.html
vector<string> tokenize(const string& str,const string& delimiters)
{
	vector<string> tokens;	
	string::size_type lastPos = 0, pos = 0;  
  	int count = 0;
  	if(str.length()<1)  return tokens;
   
  	lastPos = str.find_first_not_of(delimiters, 0); // skip delimiters at beginning.     
  	if((str.substr(0, lastPos-pos).length()) > 0)
  	{	  	
  		count = str.substr(0, lastPos-pos).length();  	
  		for(int i=0; i < count; i++)  	
  	 		tokens.push_back("");
  		if(string::npos == lastPos)
	  		tokens.push_back("");
	}

  	pos = str.find_first_of(delimiters, lastPos); // find first \"non-delimiter\".
  	while (string::npos != pos || string::npos != lastPos)
 	 {  	      	    
     		tokens.push_back( str.substr(lastPos, pos - lastPos)); // found a token, add it to the vector.				
     		lastPos = str.find_first_not_of(delimiters, pos); // skip delimiters.  Note the \"not_of\"   	   	    
		
		if((string::npos != pos) && (str.substr(pos, lastPos-pos).length() > 1))  		
  		{
  			count = str.substr(pos, lastPos-pos).length();
  			for(int i=0; i < count-1; i++)
  	 			tokens.push_back("");
		}	
  		pos = str.find_first_of(delimiters, lastPos);
	}
	return tokens;
}

void parse_sam_line(const vector<string>& tokens, vector<int>& seq_b, double& a_score, char gap_quality, int& indels, int& al_start)
{  
	seq_b.clear();
	indels = 0;
	
  	bool foundAscore = false;
  	for(int j = 11; j< tokens.size();j++)
	{
  		vector<string> tokens_as = tokenize( tokens[j].c_str(),":");
    		if( tokens_as[0] == "AS")
		{
      			a_score = atof(tokens_as[2].c_str());
      			foundAscore = true;
      			break;
    		}
 	} 
	if (!foundAscore)
	{
		cout << foundAscore << " no a_score "<< a_score << endl;
//		a_score = min_align_score_fraction * seq_b.size() *2; 
 	}
	al_start = atoi( tokens[3].c_str()) ;
 
//	if(tokens[5] != "*"){
	vector<int> isAlpha(tokens[5].size(),1);
	for(int i =0; i< tokens[5].size();i++)
	{
    		if(!isalpha(tokens[5][i]))
      			isAlpha[i] = 0;
	}
//  	}
  	vector<int> sub_length_vec;
  	vector<char> symbols;
  	int sub_length =0;
  	for(int i =0; i< tokens[5].size();i++)
	{
    		if(isAlpha[i] == 1)
		{
      			sub_length_vec.push_back( atoi(tokens[5].substr(i-sub_length,sub_length).c_str()));
      			symbols.push_back(tokens[5][i]);
      			sub_length =0;
    		}
    		if(isAlpha[i] == 0)
     			sub_length++;	
	}
  
  	int c =0;
  	for(int i =0; i< sub_length_vec.size();i++)
	{  
    		if(symbols[i] == 'S')
      			for(int j =0; j< sub_length_vec[i];j++)
				c++;
    		else if(symbols[i] == 'M')
			for(int j =0; j< sub_length_vec[i];j++)
			{
				int k = 0;
				if(tokens[9][c] == 'A' || tokens[9][c] == 'a')
			  		k=1;
				else if(tokens[9][c] == 'C' || tokens[9][c] == 'c')
			  		k=2;
				else if(tokens[9][c] == 'G' || tokens[9][c] == 'g')
			  		k=3;
				else if(tokens[9][c] == 'T' || tokens[9][c] == 't')
			  		k=4;
				seq_b.push_back(k);
				c++;
		      }
		else if(symbols[i] == 'I')
      			for(int j =0; j< sub_length_vec[i];j++)
			{	
				c++;
      			}
    		else if(symbols[i] == 'D')
      			for(int j =0; j< sub_length_vec[i];j++)
			{
				seq_b.push_back(0);
				indels++;	
      			}
	}  
}

int parseSAMpaired( string al, double max_gap_fraction, double min_align_score_fraction, double min_qual, int min_length, int max_insertln, char gap_quality, double& mean_length, vector<vector<int> >& Read_matrix, int reconstruction_start, int reconstruction_end, int& total_count, int gene_length)
{
	string line;
  	vector<string> tokens, tokens_1, tokens_2;
	int pair_counter = 0, seq_counter =0, singleton_counter = 0;
  	vector<int> seq_b, seq_b_pairs;
  	double a_score,a_score_pairs;
	int indels, indels_pairs;
  	int al_start, al_start_pairs;
    	string sRC_1, sRC_2, pairs_1, pairs_2;
  	bool part_1 = false, is_pair = false;
  	string id, id_1;
  	int RC;
    
  	ifstream inf6(al.c_str(),ios::in);
  	while( getline(inf6,line,'\n') )
	{
    		if(line[0] == '@')
      			continue;
	    	tokens = tokenize(line,"\t");
		if(tokens.size() <5)
		{
			cout << "Problem with sam file..." << endl;
			return 1;
		}
		id =  tokens[0];
	    	total_count++;
				
	    	RC =  atoi( tokens[1].c_str());
	    	stringstream strs;
	    	string sRC = binary(RC,  strs);
		int sz = sRC.size();
		if(sz > 8) // bit9-12 should be 0
	      		continue;
	    	if(sRC[sz-3] == '1' ||  sRC[sz-4]  == '1'  ) // bit 3-4 should be 0
	      		continue;
	    	if(part_1 && id == id_1 && sz == 8 && sRC[0] == '1')
		{
	      		is_pair = true;
	      		part_1 = false;
	      		pairs_2 = line;
	      		sRC_2 = sRC;
	      		tokens_2 = tokens;
	      		pair_counter++;
	   	}
	    	else
		{
	     		part_1 = true;
	      		is_pair = false;
	      		pairs_1 = line;
	      		sRC_1 = sRC;
	      		id_1 = id;
	      		tokens_1 = tokens;
	      		singleton_counter++;
	   	}
		
	    	if(is_pair)
		{
	      		parse_sam_line(tokens_1, seq_b,  a_score,  gap_quality, indels, al_start );
	      		parse_sam_line(tokens_2, seq_b_pairs,  a_score_pairs,  gap_quality, indels_pairs, al_start_pairs );
		
	      		int qual = atoi(tokens_1[4].c_str());
	      		int qual_pairs = atoi(tokens_2[4].c_str());
	
	      		if( seq_b.size() >= min_length && seq_b_pairs.size() >= min_length 
				&& qual >= min_qual && qual_pairs >= min_qual
		  		&& double(indels)/seq_b.size() < max_gap_fraction && double(indels_pairs)/seq_b_pairs.size() < max_gap_fraction 
		  		&& a_score/seq_b.size()  >  min_align_score_fraction &&  a_score_pairs/seq_b_pairs.size()  >  min_align_score_fraction)
			{	
				int StartPos = al_start;
				vector<int> SEQ_combined = seq_b;
				bool is_gap = false;
				int Nlength = 0;
	
				if(al_start_pairs < al_start)
				{
		  			StartPos = al_start_pairs;	  
		  			is_gap = false;
		  			if(al_start-(al_start_pairs +   (int)seq_b_pairs.size())>0)
		    				is_gap = true;
	  				if(is_gap)
					{
						Nlength = al_start-(al_start_pairs +   (int)seq_b_pairs.size());
	    					vector<int> Ns(Nlength,0);
						seq_b_pairs.insert(seq_b_pairs.end(), Ns.begin(), Ns.end());
						seq_b_pairs.insert(seq_b_pairs.end(), seq_b.begin(), seq_b.end());
						SEQ_combined = seq_b_pairs; 
	  				}
					else
					{
						vector<int> first_part (seq_b_pairs.begin(),seq_b_pairs.begin() + (al_start - al_start_pairs));
						first_part.insert(first_part.end(), seq_b.begin(), seq_b.end());
						SEQ_combined = first_part;
	  				}
				}
				else
				{
					if(al_start_pairs > al_start)
					{
				    		is_gap = false;
				    		if(al_start_pairs - (al_start + (int)seq_b.size()) >0)
				      			is_gap = true;
	    					if(is_gap)
						{
	      						Nlength =al_start_pairs-(al_start + (int)seq_b.size());
	      						vector<int> Ns(Nlength,0);
	      						seq_b.insert(seq_b.end(), Ns.begin(), Ns.end());
	      						seq_b.insert(seq_b.end(), seq_b_pairs.begin(), seq_b_pairs.end());
	      						SEQ_combined = seq_b;
	   		 			}
						else
						{
	      						vector<int> first_part (seq_b.begin(),seq_b.begin() + (al_start_pairs - al_start));
	      						first_part.insert(first_part.end(), seq_b_pairs.begin(), seq_b_pairs.end());
	      						SEQ_combined = first_part;
	   					 }
	  				}
				}
				if(!is_gap || Nlength < max_insertln)  
				{
					int EndPos = StartPos + SEQ_combined.size()-1; //range = reconstruction_end-reconstruction_start+1;
					if ( StartPos <= reconstruction_end && EndPos >= reconstruction_start)
					{
						vector<int> SEQ_range;
						if (StartPos < reconstruction_start)
						{
							if (EndPos <= reconstruction_end)
							{
								vector<int> SEQ_inrange(SEQ_combined.begin()+(reconstruction_start-StartPos),SEQ_combined.end());
								vector<int> Ns(gene_length-SEQ_inrange.size(),0);
								SEQ_inrange.insert(SEQ_inrange.end(),Ns.begin(),Ns.end());
								SEQ_range = SEQ_inrange;
							}
							else
							{
								vector<int> SEQ_inrange(SEQ_combined.begin()+(reconstruction_start-StartPos),SEQ_combined.begin()+(reconstruction_end-StartPos+1));
								SEQ_range = SEQ_inrange;
							}
						}
						else
						{
							if (EndPos <= reconstruction_end)
							{
								vector<int> SEQ_inrange = SEQ_combined;
								vector<int> Ns1(StartPos-reconstruction_start,0);
								SEQ_inrange.insert(SEQ_inrange.begin(),Ns1.begin(),Ns1.end());
								vector<int> Ns2(reconstruction_end-EndPos,0);
								SEQ_inrange.insert(SEQ_inrange.end(),Ns2.begin(),Ns2.end());	
								SEQ_range = SEQ_inrange;						
							}
							else
							{
								vector<int> SEQ_inrange(SEQ_combined.begin(),SEQ_combined.begin()+(reconstruction_end-StartPos+1));
								vector<int> Ns(StartPos-reconstruction_start,0);
								SEQ_inrange.insert(SEQ_inrange.begin(),Ns.begin(),Ns.end());
								SEQ_range = SEQ_inrange;
							}
						}
						Read_matrix.push_back(SEQ_range); // nReads by genome_length
			  			mean_length += SEQ_combined.size();
						seq_counter++;
					}				  
				}
	      		}
	   	 }
	}
  	return 0;
}

int parseSAM(string al, double max_gap_fraction, double min_align_score_fraction, double min_qual, int min_length, int max_insertln, char gap_quality, double& mean_length, vector<vector<int> >& Read_matrix, int reconstruction_start, int reconstruction_end, int& total_count, int gene_length)
{
	string line;
        vector<string> tokens;
	int seq_counter =0;
        vector<int> seq_b;
        double a_score;
        int indels;
        int al_start;
        string sRC_1, pairs_1;
        string id;
        int RC;
   
        ifstream inf6(al.c_str(),ios::in);
	while( getline(inf6,line,'\n') )
        {
                if(line[0] == '@')
                        continue;
                tokens = tokenize(line,"\t");
                if(tokens.size() <5)
                {
                        cout << "Problem with sam file..." << endl;
                        return 1;
                }
                id =  tokens[0];
                total_count++;
		
		RC =  atoi( tokens[1].c_str());
                stringstream strs;
                string sRC = binary(RC,  strs);
                int sz = sRC.size();
                if(sz > 8) // bit9-12 should be 0
                        continue;
		if(sRC[sz-3] == '1' ||  sRC[sz-4]  == '1'  ) // bit 3-4 should be 0
                        continue;
	//	if(sRC[sz-5] == 0 ||  sRC[sz-5]  == 1  ) 
	//	{
		parse_sam_line(tokens, seq_b,  a_score,  gap_quality, indels, al_start );
		
		int qual =  atoi( tokens[4].c_str());
		if( seq_b.size() >= min_length && qual >= min_qual && double(indels)/seq_b.size() < max_gap_fraction && a_score/seq_b.size()  >  min_align_score_fraction)
		{
			int StartPos = al_start;
			int EndPos = StartPos + seq_b.size()-1; //range = reconstruction_end-reconstruction_start+1;
			
			if ( StartPos <= reconstruction_end && EndPos >= reconstruction_start)
                        {
				vector<int> SEQ_range;
				if (StartPos < reconstruction_start)
				{
					if (EndPos <= reconstruction_end)
					{
						vector<int> SEQ_inrange(seq_b.begin()+(reconstruction_start-StartPos),seq_b.end());
                                                vector<int> Ns(gene_length-SEQ_inrange.size(),0);
                                                SEQ_inrange.insert(SEQ_inrange.end(),Ns.begin(),Ns.end());
                                                SEQ_range = SEQ_inrange;
					}
					else
					{
						vector<int> SEQ_inrange(seq_b.begin()+(reconstruction_start-StartPos),seq_b.begin()+(reconstruction_end-StartPos+1));
                                                SEQ_range = SEQ_inrange;
					}
				}
				else
				{
					if (EndPos <= reconstruction_end)
					{
						vector<int> SEQ_inrange = seq_b;
                                                vector<int> Ns1(StartPos-reconstruction_start,0);
                                                SEQ_inrange.insert(SEQ_inrange.begin(),Ns1.begin(),Ns1.end());
                                                vector<int> Ns2(reconstruction_end-EndPos,0);
                                                SEQ_inrange.insert(SEQ_inrange.end(),Ns2.begin(),Ns2.end());
                                                SEQ_range = SEQ_inrange;
                                        }
                                        else
                                        {
						vector<int> SEQ_inrange(seq_b.begin(),seq_b.begin()+(reconstruction_end-StartPos+1));
                                                vector<int> Ns(StartPos-reconstruction_start,0);
                                                SEQ_inrange.insert(SEQ_inrange.begin(),Ns.begin(),Ns.end());
                                                SEQ_range = SEQ_inrange;
                                        }
				}
			
				Read_matrix.push_back(SEQ_range); // nReads by genome_length
				mean_length += seq_b.size();
                	        seq_counter++;
			}
		}
	//	}

	}
	return 0;
}

void callSNV(vector<vector<int> >& Read_matrix, vector<vector<int> >& Allele_freq, vector<vector<int> >& SNV_matrix,  vector<vector<int> >& SNV_freq, vector<int>& Homo_seq, vector<int>& SNV_pos, int nReads, int gene_length, double SNV_thres, int& nSNV, vector<int>& deleted_reads_list)
{
	for (int j=0; j<gene_length; j++)
	{
		vector<int> count(5,0);
		for (int i=0; i<nReads; i++)
		{
			switch (Read_matrix[i][j])
			{
				case 1:
					count[0] = count[0] + 1;
					break;
				case 2:
					count[1] = count[1] + 1;
					break;
				case 3:
					count[2] = count[2] + 1;
					break;
	 			case 4:
					count[3] = count[3] + 1;
					break;			
			}
		}
		count[4] = count[0] + count[1] + count[2] +count[3];
		
		Allele_freq.push_back(count); // gene_length by 5
	}
	
	// get SNV position and corresponding allele_freq and SNV_matrix
	vector<vector<int> > temp_snv_matrix;
	for (int i = 0; i< gene_length; i++)
	{
		vector<int> allele_delete;
		int allele_sum = 0;
		int count = 0;
		if (Allele_freq[i][4] > 0)
		{
			for (int j = 0; j<4 ;j++)
			{
				double ratio = Allele_freq[i][j]/double(Allele_freq[i][4]);
				if (ratio < SNV_thres)
					allele_delete.push_back(0);
				else
				{
					allele_delete.push_back(Allele_freq[i][j]);
					allele_sum = allele_sum + Allele_freq[i][j];
					count++;
				}		
			}
			allele_delete.push_back(allele_sum);
		}
		// delete potential sequencing error
		vector<int> error_delete;
		if (count > 1)
		{
			SNV_freq.push_back(allele_delete); // nSNV by 5
			SNV_pos.push_back(i); // nSNV by 1
			
			for (int j = 0; j< nReads; j++)
			{
				if (Read_matrix[j][i]!=0)
				{
					if (allele_delete[Read_matrix[j][i]-1]!=0)
						error_delete.push_back(Read_matrix[j][i]);
					else
						error_delete.push_back(0);
				}
				else
					error_delete.push_back(0);
			}
			temp_snv_matrix.push_back(error_delete); // nSNV by nReads after error correction
			Homo_seq.push_back(0);
		}	
		else
		{
			int maxind = -1;
			int maxval = 0;
			for (int j = 0;j<4;j++)
			{
				if (maxval < Allele_freq[i][j])
				{
					maxval = Allele_freq[i][j];
					maxind = j;
				}
			}	
			Homo_seq.push_back(maxind+1);
		}
	}
	nSNV = SNV_pos.size();
	
	// delete empty fragment
	for (int i=0; i< nReads; i++)
	{ 
		int count = 0;
		vector<int> frag_delete;
		for(int j=0; j<nSNV; j++)
		{
			if (temp_snv_matrix[j][i]!=0)
				count++;
		}
		if (count!=0)
		{	
			for(int j=0; j<nSNV; j++)
				frag_delete.push_back(temp_snv_matrix[j][i]);
			SNV_matrix.push_back(frag_delete);
		}
		else
			deleted_reads_list.push_back(i);
	}
}

// initial similarity of each pair of clusters
void initial_simtable(vector<vector<int> >& simtable_pair, vector<double>& simtable_score, vector<vector<int> > SNV_matrix, int nFrag, int nSNV, double over_thres, vector <double> coverage_rate)
{	
	for (int i=0; i< nFrag; i++)
	{
		simtable_pair[i][0] = i;
		simtable_pair[i][1] = -1;
	}
	
	for (int i = 0; i < nFrag; i++)
	{
		if (i%10000 == 0)
			cout << i << endl;
		double max_score = -100;
		int max_pair = -1;
		for (int j=0; j<nFrag; j++)
		{
			if (i!=j)
			{
				int L_min=0, L=0, Lo=0, La=0, Lb=0;
				int conf_index = 0; // 0: no conflict, 1: conflict
				vector<int> pos_index(nSNV,0);
				for (int n=0; n<nSNV; n++)
				{
					if (SNV_matrix[i][n]>0 && SNV_matrix[j][n]>0)
					{
						L++, Lo++, La++, Lb++;
						pos_index[n] = 3;
						if (SNV_matrix[i][n]!=SNV_matrix[j][n])
						{
							conf_index = 1;
							break;		
						}
					}
					else if (SNV_matrix[i][n]>0 && SNV_matrix[j][n]==0) 
					{
						pos_index[n] = 1;
						L++, La++;
					}
					else if (SNV_matrix[i][n]==0 && SNV_matrix[j][n]>0) 
					{
						pos_index[n] = 2;
						L++, Lb++;
					}
				}
				L_min = min(La,Lb);
				if (conf_index==0 && Lo/double(L_min)>=0.5)
				{
					if (Lo/double(L)>=over_thres || Lo==L_min)
					{
						vector<double> score(L,0);
						int score_index=-1;
						for (int n=0; n<nSNV; n++)
						{
							if (pos_index[n]==3)
							{
								score_index++;
								score[score_index]=1;
							}
							else if (pos_index[n]==1 || pos_index[n]==2)
							{
								score_index++;
								score[score_index] = (0.5-coverage_rate[n])/(1-coverage_rate[n]);
							}
						}
						double sim_score = accumulate(score.begin(),score.end(),0)/double(L); 
						if (sim_score>max_score)
						{
							max_score=sim_score;
							max_pair=j;
							
							if (max_score==1)
								break;	
						}
					}
				}
			}
		}
		simtable_score[i] = max_score;
		simtable_pair[i][1] = max_pair;
	}
}

int find_maxscore(vector<double> simtable_score, int nNode)
{
	int index = -1;
	double max_score = -100;
	for (int i=0; i<nNode; i++)
	{
		if (simtable_score[i] == 1)
		{
			index = i;
			break;
		}
		else
		{
			if (simtable_score[i]>max_score)
			{
				index = i;
				max_score = simtable_score[i];
			}
		}
	}
	return index;
}

int find_pair(int b, vector<vector<int> > simtable_pair, int nNode)
{
	int pair_row;
	for (int i=0; i<nNode; i++)
	{
		if (simtable_pair[i][0]==b)
			pair_row = i;
	}
	return pair_row;
}

int find_group(int a, int nGroup, vector<vector<int> > group_matrix)
{
	int grpind;
	for (int n=0; n<nGroup; n++)
	{
		if (a==group_matrix[n][0])
		{
			grpind = n;
			break;
		}
	}
	return grpind;
}

// update the measure of similarity after merging two clusters 
// use partial maximum array algorithm
void update_simtable(vector<vector<int> >& simtable_pair, vector<double>& simtable_score, vector<vector<int> > SNV_matrix, int nNode, int nSNV, double over_thres, vector <double> coverage_rate, vector<vector<int> > group_matrix, vector<vector<int> > superread_matrix, vector<vector<int> > count_matrix, vector<vector<double> > score_matrix, vector<int> nMember, int nGroup, vector<int> assigned_index, int grp_ind, int ab, int a, int b, int sim_row)
{
	vector<double> temp_max_score(nNode,-100);

	double max_score = -100;
	int max_pair = -1;
	
	vector<int> seq1;
	seq1.assign(superread_matrix[grp_ind].begin(), superread_matrix[grp_ind].end());
	vector<int> count1;
	count1.assign(count_matrix[grp_ind].begin(), count_matrix[grp_ind].end());
	for (int j=0; j<nNode; j++)
	{
		if (simtable_pair[j][0]!=-1 && simtable_pair[j][0]!=ab)
		{
			int r = simtable_pair[j][0];
			int L_min=0, L=0, Lo=0, La=0, Lb=0;
			int conf_index = 0; // 0: no conflict, 1: conflict
			vector<int> pos_index(nSNV,0);
			
			vector<int> seq2;
			vector<int> count2;
			int ntotalMember;
			if (assigned_index[r] == 1)
			{
				int grp_row;
				for (int n=0; n<nGroup; n++)
				{
					if (r==group_matrix[n][0])
					{
						grp_row = n;
						break;
					}
				}
				seq2.assign(superread_matrix[grp_row].begin(), superread_matrix[grp_row].end());
				count2.assign(count_matrix[grp_row].begin(), count_matrix[grp_row].end());
				ntotalMember = nMember[grp_ind] + nMember[grp_row];
			}
			else
			{
				seq2.assign(SNV_matrix[r].begin(),SNV_matrix[r].end());	
				for (int n=0; n<seq2.size(); n++)
				{
					if (seq2[n]==0)
						count2.push_back(0);
					else
						count2.push_back(1);
				}
				ntotalMember = nMember[grp_ind] + 1;
			}
			for (int n = 0; n<nSNV; n++)
			{	
				if (seq1[n]>0 && seq2[n]>0)
				{
					L++, Lo++, La++, Lb++;
					pos_index[n] = 3;
					if (seq1[n]!=seq2[n])
					{
						conf_index = 1;
						break;	
					}
				}
				else if (seq1[n]>0 && seq2[n]==0)
				{
					pos_index[n] = 1;
					L++, La++;
				}
				else if (seq1[n]==0 && seq2[n]>0)
				{
					pos_index[n] = 2;
					L++, Lb++;
				}
			}
			L_min = min(La,Lb);
			if (conf_index==0 && Lo/double(L_min)>=0.5)
			{
				if (Lo/double(L)>=over_thres || Lo==L_min)
				{
					vector<double> score(L,0);
					int score_index=-1;
					for (int n=0; n<nSNV; n++)
					{
						if (pos_index[n]!=0)
						{
							score_index++;
							double grp_cov = (count1[n]+count2[n])/double(ntotalMember);
							score[score_index] = (grp_cov-coverage_rate[n])/(1-coverage_rate[n]);
						}
					}
					double sim_score = accumulate(score.begin(),score.end(),0)/double(L); 
					temp_max_score[j] = sim_score; 
					if (sim_score>max_score)
					{
						max_score=sim_score;
						max_pair=r;
					}
				}
			}
		}
	}
	simtable_score[sim_row] = max_score;
	simtable_pair[sim_row][1] = max_pair;
	
	for (int j=0; j<nNode; j++)
	{	
		if (simtable_pair[j][0]!=-1 && simtable_pair[j][0]!=ab)
		{
			int a1 = simtable_pair[j][0];
			if (temp_max_score[j] > simtable_score[j])
			{
				simtable_score[j] = temp_max_score[j];
				simtable_pair[j][1] = ab;
			}
			else if (temp_max_score[j] == simtable_score[j] && temp_max_score[j] > -100)
			{
				simtable_pair[j][1] = ab;
			}
			else if (temp_max_score[j] < simtable_score[j])
			{
				if (simtable_pair[j][1] == a || simtable_pair[j][1] == b)
				{
					vector<double> temp_max_score1(nNode,-100);
					double max_score1 = -100;
					int max_pair1 = -1;
				
					vector<int> seq1;
					vector<int> count1;
					int ntotalMember1;
					if (assigned_index[a1] == 1)
					{
						int grp_row;
						for (int n=0; n<nGroup; n++)
						{
							if (a1==group_matrix[n][0])
							{
								grp_row = n;
								break;
							}
						}
						seq1.assign(superread_matrix[grp_row].begin(), superread_matrix[grp_row].end());
						count1.assign(count_matrix[grp_row].begin(), count_matrix[grp_row].end());
						ntotalMember1= nMember[grp_row];
					}
					else
					{
						seq1.assign(SNV_matrix[a1].begin(),SNV_matrix[a1].end());	
						for (int n=0; n<seq1.size(); n++)
						{
							if (seq1[n]==0)
								count1.push_back(0);
							else
								count1.push_back(1);
						}
						ntotalMember1 = 1;
					}
					for (int j1=0; j1<nNode; j1++)
					{
						if (simtable_pair[j1][0]!=-1 && simtable_pair[j1][0]!=a1)
						{
							int b1 = simtable_pair[j1][0];
							
							
							vector<int> seq2;
							vector<int> count2;
							int ntotalMember;
							if (assigned_index[b1] == 1)
							{
								int grp_row;
								for (int n=0; n<nGroup; n++)
								{
									if (b1==group_matrix[n][0])
									{
										grp_row = n;
										break;
									}
								}
								seq2.assign(superread_matrix[grp_row].begin(), superread_matrix[grp_row].end());
								count2.assign(count_matrix[grp_row].begin(), count_matrix[grp_row].end());
								ntotalMember = ntotalMember1 + nMember[grp_row];
							}
							else
							{
								seq2.assign(SNV_matrix[b1].begin(),SNV_matrix[b1].end());	
								for (int n=0; n<seq2.size(); n++)
								{
									if (seq2[n]==0)
										count2.push_back(0);
									else
										count2.push_back(1);
								}
								ntotalMember = ntotalMember1 + 1;
							}
							int L_min=0, L=0, Lo=0, La=0, Lb=0;
							int conf_index = 0; // 0: no conflict, 1: conflict
							vector<int> pos_index(nSNV,0);
							for (int n = 0; n<nSNV; n++)
							{	
								if (seq1[n]>0 && seq2[n]>0)
								{
									L++, Lo++, La++, Lb++;
									pos_index[n] = 3;
									if (seq1[n]!=seq2[n])
									{
										conf_index = 1;
										break;	
									}
								}
								else if (seq1[n]>0 && seq2[n]==0)
								{
									pos_index[n] = 1;
									L++, La++;
								}
								else if (seq1[n]==0 && seq2[n]>0)
								{
									pos_index[n] = 2;
									L++, Lb++;
								}
							}
							L_min = min(La,Lb);
							if (conf_index==0 && Lo/double(L_min)>=0.5)
							{
								if (Lo/double(L)>=over_thres || Lo==L_min)
								{
									vector<double> score(L,0);
									int score_index=-1;
									for (int n=0; n<nSNV; n++)
									{
										if (pos_index[n]!=0)
										{
											score_index++;
											double grp_cov = (count1[n]+count2[n])/double(ntotalMember);
											score[score_index] = (grp_cov-coverage_rate[n])/(1-coverage_rate[n]);
										}
									}
									double sim_score1= accumulate(score.begin(),score.end(),0)/double(L); 
									
									temp_max_score1[j1] = sim_score1; //////////////////// omitted/////////
									if (sim_score1>max_score1)
									{
										max_score1=sim_score1;
										max_pair1=b1;
										
										if (max_score1 == 1)
											break;
									}
								}
							}
						}
					}
					simtable_score[j] = max_score1;
					simtable_pair[j][1] = max_pair1;
				}
			}
		}
	}
}

// merge a pair of clusters and update the new cluster
void merge_group(vector<vector<int> >& simtable_pair, vector<double>& simtable_score, vector<vector<int> > SNV_matrix, int nNode, int nSNV, double over_thres, vector <double> coverage_rate, vector<vector<int> >& group_matrix, vector<vector<int> >& superread_matrix, vector<vector<int> >& count_matrix, vector<vector<double> >& score_matrix, vector<int>& nMember, int& nGroup, vector<int>& assigned_index, vector<vector<double> >& overlap_ratio, int& grp_ind, int& ab, int& a, int& b, int& sim_row, int& iter)
{
	for (iter=0;iter<nNode;iter++)
	{

		if (iter%1000 == 0)
			cout << iter << endl;
			
		int Ind = find_maxscore(simtable_score, nNode);
		if (Ind == -1)
			break;

		double Val = simtable_score[Ind];  

		a = simtable_pair[Ind][0];
		b = simtable_pair[Ind][1];
				
		int row1 = Ind;
		int row2 = find_pair(b, simtable_pair, nNode);
		
		if (simtable_score[row1]!=simtable_score[row2])
		{
			cout << "Problem with similarity measure.." << endl;
			break;
		}
		
		if (assigned_index[a]==0 && assigned_index[b]==0)
		{
			grp_ind = nGroup; 
			nGroup++;
			ab = a;
			vector<int> member(2);
			member[0] = a;
			member[1] = b;
			group_matrix.push_back(member);
			nMember.push_back(2);
			sim_row = row1;
			simtable_score[row2] = -100;
			simtable_pair[row2][0] = -1;
			simtable_pair[row2][1] = -1;
			
			int L_min=0, L=0, Lo=0, La=0, Lb=0;
			vector<int> superread_vec(nSNV,0);
			vector<int> count_vec(nSNV,0);
			vector<double> score_vec(nSNV,-100);
			for (int n=0; n<nSNV; n++)
			{
				if (SNV_matrix[a][n]>0 && SNV_matrix[b][n]>0)
				{
					L++, Lo++, La++, Lb++;
					superread_vec[n] = SNV_matrix[a][n];
					count_vec[n] = 2;
					score_vec[n] = 1;
				}
				else if (SNV_matrix[a][n]>0 && SNV_matrix[b][n]==0) 
				{
					L++, La++;
					superread_vec[n] = SNV_matrix[a][n];
					count_vec[n] = 1;
					score_vec[n] = (0.5-coverage_rate[n])/(1-coverage_rate[n]);
				}
				else if (SNV_matrix[a][n]==0 && SNV_matrix[b][n]>0) 
				{
					L++, Lb++;
					superread_vec[n] = SNV_matrix[a][n];
					count_vec[n] = 1;
					score_vec[n] = (0.5-coverage_rate[n])/(1-coverage_rate[n]);
				}
			}
			superread_matrix.push_back(superread_vec);
			count_matrix.push_back(count_vec);
			score_matrix.push_back(score_vec);
			assigned_index[a] = 1;
			assigned_index[b] = 1;
				
				
			L_min = min(La,Lb);		
			vector<double> overlap_vec(3);
			overlap_vec[0] = Lo/double(L);
			overlap_vec[1] = Lo/double(L_min);
			overlap_vec[2] = Val;
			overlap_ratio.push_back(overlap_vec);
		}
		else if (assigned_index[a]==1 && assigned_index[b]==0)
		{
			grp_ind = find_group(a, nGroup, group_matrix);
			ab = a;
			nMember[grp_ind]++;
			
			group_matrix[grp_ind].resize(nMember[grp_ind]);
			group_matrix[grp_ind][nMember[grp_ind]-1] = b;

			sim_row = row1;
			simtable_score[row2] = -100;
			simtable_pair[row2][0] = -1;
			simtable_pair[row2][1] = -1;
			
			int L_min=0, L=0, Lo=0, La=0, Lb=0;
			for (int n=0; n<nSNV; n++)
			{
				if (superread_matrix[grp_ind][n]>0 && SNV_matrix[b][n]>0)
				{
					L++, Lo++, La++, Lb++;
					count_matrix[grp_ind][n]++;
					score_matrix[grp_ind][n] = (count_matrix[grp_ind][n]/double(nMember[grp_ind])-coverage_rate[n])/(1-coverage_rate[n]);
				}
				else if (superread_matrix[grp_ind][n]>0 && SNV_matrix[b][n]==0) 
				{
					L++, La++;
					score_matrix[grp_ind][n] = (count_matrix[grp_ind][n]/double(nMember[grp_ind])-coverage_rate[n])/(1-coverage_rate[n]);
				}
				else if (superread_matrix[grp_ind][n]==0 && SNV_matrix[b][n]>0) 
				{
					L++, Lb++;
					superread_matrix[grp_ind][n] = SNV_matrix[b][n];
					count_matrix[grp_ind][n] = 1;
					score_matrix[grp_ind][n] = (1/double(nMember[grp_ind])-coverage_rate[n])/(1-coverage_rate[n]);
				}
			}
			assigned_index[b] = 1;
			
			L_min = min(La,Lb);		
			vector<double> overlap_vec(3);
			overlap_vec[0] = Lo/double(L);
			overlap_vec[1] = Lo/double(L_min);
			overlap_vec[2] = Val;
			overlap_ratio.push_back(overlap_vec);
		}
		else if (assigned_index[a]==0 && assigned_index[b]==1)
		{
		
			grp_ind = find_group(b, nGroup, group_matrix);
			ab = b;
			nMember[grp_ind]++;
			
			group_matrix[grp_ind].resize(nMember[grp_ind]);
			group_matrix[grp_ind][nMember[grp_ind]-1] = a;

			sim_row = row2;
			simtable_score[row1] = -100;
			simtable_pair[row1][0] = -1;
			simtable_pair[row1][1] = -1;
			
			int L_min=0, L=0, Lo=0, La=0, Lb=0;
			for (int n=0; n<nSNV; n++)
			{
				if (SNV_matrix[a][n]>0 && superread_matrix[grp_ind][n]>0)
				{
					L++, Lo++, La++, Lb++;
					count_matrix[grp_ind][n]++;
					score_matrix[grp_ind][n] = (count_matrix[grp_ind][n]/double(nMember[grp_ind])-coverage_rate[n])/(1-coverage_rate[n]);
				}
				else if (SNV_matrix[a][n]>0 && superread_matrix[grp_ind][n]==0) 
				{
					L++, La++;
					superread_matrix[grp_ind][n] = SNV_matrix[a][n];
					count_matrix[grp_ind][n] = 1;
					score_matrix[grp_ind][n] = (1/double(nMember[grp_ind])-coverage_rate[n])/(1-coverage_rate[n]);
				}
				else if (SNV_matrix[a][n]==0 && superread_matrix[grp_ind][n]>0) 
				{
					L++, Lb++;
					score_matrix[grp_ind][n] = (count_matrix[grp_ind][n]/double(nMember[grp_ind])-coverage_rate[n])/(1-coverage_rate[n]);
				}
			}
			assigned_index[a] = 1;
					
			L_min = min(La,Lb);		
			vector<double> overlap_vec(3);
			overlap_vec[0] = Lo/double(L);
			overlap_vec[1] = Lo/double(L_min);
			overlap_vec[2] = Val;
			overlap_ratio.push_back(overlap_vec);
		}
		else if (assigned_index[a]==1 && assigned_index[b]==1)
		{
			int grp_inda = find_group(a, nGroup, group_matrix);
			int grp_indb = find_group(b, nGroup, group_matrix);
			if (grp_inda > grp_indb)
			{
				swap(a,b);
				swap(grp_inda,grp_indb);
				swap(row1,row2);				
			}
			grp_ind = grp_inda;
			
			group_matrix[grp_inda].insert(group_matrix[grp_inda].end(),group_matrix[grp_indb].begin(),group_matrix[grp_indb].end());
			group_matrix[grp_indb].assign(group_matrix[nGroup-1].begin(),group_matrix[nGroup-1].end());
			group_matrix.erase(group_matrix.begin()+nGroup-1);
			
			ab = group_matrix[grp_inda][0]; // ab=a;
			
			nMember[grp_inda] = group_matrix[grp_inda].size();
			nMember[grp_indb] = group_matrix[grp_indb].size();
			
			sim_row = row1;
			simtable_score[row2] = -100;
			simtable_pair[row2][0] = -1;
			simtable_pair[row2][1] = -1;
			
			int L_min=0, L=0, Lo=0, La=0, Lb=0;
			for (int n=0; n<nSNV; n++)
			{
				if (superread_matrix[grp_inda][n]>0 && superread_matrix[grp_indb][n]>0)
				{
					L++, Lo++, La++, Lb++;
					
					count_matrix[grp_inda][n] = count_matrix[grp_inda][n] + count_matrix[grp_indb][n];
					score_matrix[grp_inda][n] = (count_matrix[grp_inda][n]/double(nMember[grp_inda])-coverage_rate[n])/(1-coverage_rate[n]);
					
					superread_matrix[grp_indb][n] = 0;   ///////////unnecessary
					count_matrix[grp_indb][n] = 0;  ///////////unnecessary
					score_matrix[grp_indb][n] = -100;  ///////////unnecessary
				}
				else if (superread_matrix[grp_inda][n]>0 && superread_matrix[grp_indb][n]==0) 
				{
					L++, La++;
					score_matrix[grp_inda][n] = (count_matrix[grp_inda][n]/double(nMember[grp_inda])-coverage_rate[n])/(1-coverage_rate[n]);
				}
				else if (superread_matrix[grp_inda][n]==0 && superread_matrix[grp_indb][n]>0) 
				{
					L++, Lb++;
					superread_matrix[grp_inda][n] = superread_matrix[grp_indb][n];
					count_matrix[grp_inda][n] = count_matrix[grp_indb][n];
					score_matrix[grp_inda][n] = (count_matrix[grp_inda][n]/double(nMember[grp_inda])-coverage_rate[n])/(1-coverage_rate[n]);
					
					superread_matrix[grp_indb][n] = 0;   ///////////unnecessary
					count_matrix[grp_indb][n] = 0;    ///////////unnecessary
					score_matrix[grp_indb][n] = -100;   ///////////unnecessary
				}
			}
			superread_matrix[grp_indb].assign(superread_matrix[nGroup-1].begin(),superread_matrix[nGroup-1].end());
			superread_matrix.erase(superread_matrix.begin()+nGroup-1);
			count_matrix[grp_indb].assign(count_matrix[nGroup-1].begin(),count_matrix[nGroup-1].end());
			count_matrix.erase(count_matrix.begin()+nGroup-1);
			score_matrix[grp_indb].assign(score_matrix[nGroup-1].begin(),score_matrix[nGroup-1].end());
			score_matrix.erase(score_matrix.begin()+nGroup-1);			
			nMember.erase(nMember.begin()+nGroup-1);
			nGroup--;
			
			L_min = min(La,Lb);		
			vector<double> overlap_vec(3);
			overlap_vec[0] = Lo/double(L);
			overlap_vec[1] = Lo/double(L_min);
			overlap_vec[2] = Val;
			overlap_ratio.push_back(overlap_vec);
		}
		
		
		update_simtable(simtable_pair, simtable_score, SNV_matrix, nNode, nSNV, over_thres, coverage_rate, group_matrix, superread_matrix, count_matrix, score_matrix, nMember, nGroup, assigned_index, grp_ind, ab, a, b, sim_row);
		
	}
}

// measure the new similarity of pair-wise clusters after adjucting overlap_ratio thresthold
void get_simtable(vector<int> nodelist, vector<vector<int> >& simtable_pair, vector<double>& simtable_score, vector<vector<int> > SNV_matrix, int nNode, int nSNV, double over_thres, vector <double> coverage_rate, vector<int> assigned_index, vector<vector<int> > group_matrix, vector<vector<int> > superread_matrix, vector<vector<int> > count_matrix, vector<int> nMember, int nGroup)
{
	for (int i=0; i<nNode; i++)
	{
		simtable_pair[i][0] = nodelist[i];
		simtable_pair[i][1] = -1;
		
		int a = nodelist[i];
		double max_score = -100;
		int max_pair = -1;
		
		vector<int> seq1;
		vector<int> count1;
		int ntotalMember1;
		if (assigned_index[a] == 1)
		{
			int grp_row;
			for (int n=0; n<nGroup; n++)
			{
				if (a==group_matrix[n][0])
				{
					grp_row = n;
					break;
				}
			}
			seq1.assign(superread_matrix[grp_row].begin(), superread_matrix[grp_row].end());
			count1.assign(count_matrix[grp_row].begin(), count_matrix[grp_row].end());
			ntotalMember1= nMember[grp_row];
		}
		else
		{
			seq1.assign(SNV_matrix[a].begin(),SNV_matrix[a].end());	
			for (int n=0; n<seq1.size(); n++)
			{
				if (seq1[n]==0)
					count1.push_back(0);
				else
					count1.push_back(1);
			}
			ntotalMember1 = 1;
		}
		for (int j=0; j<nNode; j++)
		{
			int b = nodelist[j];
			if (a!=b)
			{
				vector<int> seq2;
				vector<int> count2;
				int ntotalMember;
				if (assigned_index[b] == 1)
				{
					int grp_row;
					for (int n=0; n<nGroup; n++)
					{
						if (b==group_matrix[n][0])
						{
							grp_row = n;
							break;
						}
					}
					seq2.assign(superread_matrix[grp_row].begin(), superread_matrix[grp_row].end());
					count2.assign(count_matrix[grp_row].begin(), count_matrix[grp_row].end());
					ntotalMember = ntotalMember1 + nMember[grp_row];
				}
				else
				{
					seq2.assign(SNV_matrix[b].begin(),SNV_matrix[b].end());	
					for (int n=0; n<seq2.size(); n++)
					{
						if (seq2[n]==0)
							count2.push_back(0);
						else
							count2.push_back(1);
					}
					ntotalMember = ntotalMember1 + 1;
				}
				int L_min=0, L=0, Lo=0, La=0, Lb=0;
				int conf_index = 0; // 0: no conflict, 1: conflict
				vector<int> pos_index(nSNV,0);
				for (int n = 0; n<nSNV; n++)
				{	
					if (seq1[n]>0 && seq2[n]>0)
					{
						L++, Lo++, La++, Lb++;
						pos_index[n] = 3;
						if (seq1[n]!=seq2[n])
						{
							conf_index = 1;
							// sim_score = -100;
							break;	
						}
					}
					else if (seq1[n]>0 && seq2[n]==0)
					{
						pos_index[n] = 1;
						L++, La++;
					}
					else if (seq1[n]==0 && seq2[n]>0)
					{
						pos_index[n] = 2;
						L++, Lb++;
					}
				}	
				L_min = min(La,Lb);
				if (conf_index==0 && Lo/double(L_min)>=0.5)
				{
					if (Lo/double(L)>=over_thres || Lo==L_min)
					{
						vector<double> score(L,0);
						int score_index=-1;
						for (int n=0; n<nSNV; n++)
						{
							if (pos_index[n]!=0)
							{
								score_index++;
								double grp_cov = (count1[n]+count2[n])/double(ntotalMember);
								score[score_index] = (grp_cov-coverage_rate[n])/(1-coverage_rate[n]);
							}
						}
						double sim_score= accumulate(score.begin(),score.end(),0)/double(L); 	
						if (sim_score>max_score)
						{
							max_score=sim_score;
							max_pair=b;
							
							if (max_score==1)
								break;	
						}
					}
				}
			}
		} 
		simtable_score[i] = max_score;
		simtable_pair[i][1] = max_pair;
	}	
}

// aggromerative clustering for super-reads construction
void aggromerative_clustering(vector<vector<int> > SNV_matrix, vector<vector<int> > SNV_freq, int nFrag, int nSNV, double overlap_ratio_start, double overlap_ratio_end, string zonename, int min_length, vector<vector<int> >& group_matrix, vector<vector<int> >& superread_matrix, vector<vector<int> >& count_matrix, vector<vector<double> >& score_matrix, vector<int>& nMember, int & nGroup)
{
	vector <double> coverage_rate(nSNV); // SNV position-wise coverage rate
	for (int i=0; i<nSNV; i++)
		coverage_rate[i] = SNV_freq[i][4]/double(nFrag);
	
	cout << "Iteration start from overlap_ratio " << overlap_ratio_start << endl;
	
	clock_t start, end;
	double cpu_time_used;
	start = clock();
	vector<vector<int> > simtable_pair(nFrag,vector<int>(2));
	vector<double> simtable_score(nFrag,-100);
	double over_thres = overlap_ratio_start;
	initial_simtable(simtable_pair, simtable_score, SNV_matrix, nFrag, nSNV, over_thres, coverage_rate);
       	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for getting initial simtable: "  << cpu_time_used << endl << endl;
	cout << "nSNV: " <<nSNV<< endl; 
	
	start = clock();
	int nNode = nFrag;
	vector<int> assigned_index(nNode,0);
	vector<vector<double> > overlap_ratio;	
	int grp_ind, ab, a, b, sim_row, iter;	
	merge_group(simtable_pair, simtable_score, SNV_matrix, nNode, nSNV, over_thres, coverage_rate, group_matrix, superread_matrix, count_matrix, score_matrix, nMember, nGroup, assigned_index, overlap_ratio, grp_ind, ab, a, b, sim_row, iter);
	cout << "After "<<iter<<" iterations with overlap_ratio "<<over_thres << ", " << nGroup <<" superreads are generated."<<endl;
	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for first iteration: "  << cpu_time_used << endl << endl;

	while(over_thres>overlap_ratio_end+0.1)
	{	
		start = clock();
		over_thres = over_thres - 0.1;
		
		vector<int> nodelist;
		for (int i=0; i<nNode;i++)
		{
			if (simtable_pair[i][0]!=-1)
				nodelist.push_back(simtable_pair[i][0]);
		}
		nNode = nodelist.size();
		
		simtable_pair.clear();
		simtable_score.clear();		
		simtable_score.assign(nNode,-100);
		simtable_pair.assign(nNode, vector<int>(2,-1));
		
		cout << "Clustering for overlap_ratio " << over_thres << " starts with " << nNode<< " nodes. " <<endl;
		get_simtable(nodelist, simtable_pair, simtable_score, SNV_matrix, nNode, nSNV, over_thres, coverage_rate, assigned_index, group_matrix, superread_matrix, count_matrix, nMember, nGroup);
		merge_group(simtable_pair, simtable_score, SNV_matrix, nNode, nSNV, over_thres, coverage_rate, group_matrix, superread_matrix, count_matrix, score_matrix, nMember, nGroup, assigned_index, overlap_ratio, grp_ind, ab, a, b, sim_row, iter);
		cout << "After "<<iter<<" iterations, " << nGroup <<" superreads remain."<<endl;
		end = clock();	
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		cout << "CPU time for this iteration: "  << cpu_time_used << endl<<endl;
	}
}

// enumerate all possible extension of strains
void getHapSeq(vector<vector<int> >& cand, int& numcand, int M, vector<vector<int> > super_mat, vector<vector<int> > pre_cand, int pre_numcand, int s, vector<vector<int> > count_mat, vector<int> permu, vector<vector<double> > errtable, int K, vector<int>& preseq_index, vector<int> remain, vector<int>& vec_relNcov, vector<double>& vec_relPR)
{
	int totalcov = 0;
	int Npermu = permu.size();
	vector<vector<double> > relPR(pre_numcand,vector<double>(Npermu));
	vector<vector<int> > relNcov(pre_numcand,vector<int>(Npermu));
	vector<vector<int> > final_numpermucand(pre_numcand,vector<int>(Npermu));
	
	vector<double> PR(Npermu,1);
	int num = 0;
	for (int i=0; i<M; i++)
	{
		if (super_mat[i][s]!=0)
		{
			num++;
			for (int p=0; p<Npermu; p++)
			{
				if (permu[p] == super_mat[i][s])
				{
					PR[p] = PR[p]*errtable[i][s];
				}
				else
				{
					PR[p] = PR[p]*(1-errtable[i][s]);
				}
			}
		}
	}
	vector<double> allreads_PR(Npermu,1);
	vector<double> allreads_PR_relPR(Npermu,1);
	double invnum = 1/double(num);
	
	for (int i=0; i<Npermu; i++)
		allreads_PR[i] = pow(PR[i],invnum);
	double sum_allreads_PR = accumulate(allreads_PR.begin(),allreads_PR.end(),0.0);
	for (int i=0; i<Npermu; i++)
		allreads_PR_relPR[i] = allreads_PR[i]/double(sum_allreads_PR);

	for (int c=0; c<pre_numcand; c++)
	{
		int h_pre[s];
		for (int i=0; i<s; i++)
		{	
			h_pre[i] = pre_cand[c][i];
		}
		vector<double> PR2(Npermu,1);
		vector<int> Ncov(Npermu,0);
		int num_msr = 0;
		vector<int> num_permucand(Npermu,0);
		vector<int> matchread;
		for (int i=0; i<M; i++)
		{
			if (super_mat[i][s]!=0)
			{
				int numzero = 0;
				int unmatchind = 0;
				for (int j=0; j<s; j++)
				{
					if (super_mat[i][j]!=0)
					{
						if (super_mat[i][j]!=h_pre[j])
						{
							unmatchind = 1;
							break;
						}
					}
					else
						numzero++;
				}
				if (numzero!=s && unmatchind==0)
				{
					num_msr++;
					matchread.push_back(i);
					for (int p=0; p<Npermu; p++)
					{
						if (permu[p] == super_mat[i][s])
						{
							num_permucand[p]++;
							PR2[p] = PR2[p]*errtable[i][s];
							Ncov[p] = Ncov[p] + count_mat[i][s];
						}
						else
						{
							PR2[p] = PR2[p]*(1-errtable[i][s]);
						}
					}
				}
			}
		}
		totalcov = totalcov + accumulate(Ncov.begin(),Ncov.end(),0);
		if (num_msr == 0)
		{
			for (int i=0; i<Npermu; i++)
				relPR[c][i] = allreads_PR_relPR[i];
		}
		else
		{
			double sumPR2 = accumulate(PR2.begin(),PR2.end(),0.0);
			for (int i=0; i<Npermu; i++)
				relPR[c][i] = PR2[i]/double(sumPR2);
		}
		for (int i=0; i<Npermu; i++)
			relNcov[c][i] = Ncov[i];
		for (int i=0; i<Npermu; i++)
			final_numpermucand[c][i] = num_permucand[i];

	}
	for (int c=0; c<pre_numcand; c++)
	{
		vector<int> permucand;
		for (int i=0; i<Npermu; i++)
		{
			if (final_numpermucand[c][i]!=0)
				permucand.push_back(final_numpermucand[c][i]);
		}
		if (permucand.size()==0)
		{	
			for (int p=0; p<Npermu; p++)
			{
				for (int k=0; k<s; k++)
					cand[numcand][k] = pre_cand[c][k];
				cand[numcand][s] = permu[p];
				preseq_index[numcand] = remain[c];
				vec_relNcov.push_back(relNcov[c][p]);
				vec_relPR.push_back(relPR[c][p]);	
				numcand++;				
			}
		}
		else
		{
			for (int p=0; p<Npermu; p++)
			{
				if (relPR[c][p]>=delta0 || relNcov[c][p]>=totalcov*0.01)
				{
					for (int k=0; k<s; k++)
						cand[numcand][k] = pre_cand[c][k];
					cand[numcand][s] = permu[p];
					preseq_index[numcand] = remain[c];
					vec_relNcov.push_back(relNcov[c][p]);
					vec_relPR.push_back(relPR[c][p]);
					numcand++;
				}
			}
		}
	}
	cand.erase(cand.begin()+numcand,cand.end());
	preseq_index.erase(preseq_index.begin()+numcand,preseq_index.end());
}

void minIndVal(vector<int> vec, int& minval, int& minind)
{	
	minval = vec[0];
	minind = 0;
	for (int i=1; i<vec.size(); i++)
	{
		if (vec[i]<minval)
		{
			minval = vec[i];
			minind = i;
		}
	}
}
void minIndValdouble(vector<double> vec, double& minval, int& minind)
{	
	minval = vec[0];
	minind = 0;
	for (int i=1; i<vec.size(); i++)
	{
		if (vec[i]<minval)
		{
			minval = vec[i];
			minind = i;
		}
	}
}

void maxIndValdouble(vector<double> vec, double& maxval, int& maxind)
{	
	maxval = vec[0];
	maxind = 0;
	for (int i=1; i<vec.size(); i++)
	{
		if (vec[i]>maxval)
		{
			maxval = vec[i];
			maxind = i;
		}
	}
}

void go(int offset, int k, vector<vector<int> >& allcomb, vector<int> list, vector<int>& combnation) 
{
	if (k == 0) 
	{
		allcomb.push_back(combnation);
    		return;
  	}
  	for (int i = offset; i <= list.size() - k; ++i) 
	{
    		combnation.push_back(list[i]);
    		go(i+1, k-1, allcomb, list, combnation);
    		combnation.pop_back();
  	}
}

void combnk(int n, int k,vector<vector<int> >& comb) 
{	
	vector<vector<int> > allcomb;
	vector<int> list;
	vector<int> combnation;

  	for (int i = 0; i < n; ++i) { list.push_back(i); }

  	go(0, k, allcomb,list, combnation);
	comb = allcomb;
}

void getLogLikeli(vector<double>& logPR, int numhapcand, int M, int K, int s, vector<vector<int> > superread_mat, vector<vector<vector<int> > > hap, vector<vector<double> > errtable)
{
	for (int c=0; c<numhapcand; c++)
	{
		vector<double> Pr(M,0);
		vector<double> logPr(M,0);
		for (int i=0; i<M; i++)
		{
			vector<double> pr(K,1);
			for (int k=0; k<K; k++)
			{
				for (int j=0; j<=s; j++)
				{
					if (superread_mat[i][j]!=0)
					{
						if (superread_mat[i][j] == hap[c][k][j])
							pr[k] = pr[k]*errtable[i][j];
						else
							pr[k] = pr[k]*(1-errtable[i][j]);
					}
				}
			}
			Pr[i] = accumulate(pr.begin(),pr.end(),0.0)/double(K);
			logPr[i] = log(Pr[i]);
		}
		logPR[c] = accumulate(logPr.begin(),logPr.end(),0.0);
	}

}

void getMECrate(vector<int>& MEC, vector<double>& MECrate, int numhapcand, int M, int K, int s, vector<vector<int> > superread_mat, vector<vector<vector<int> > > hap, vector<vector<int> > count_mat)
{
	for (int c=0; c<numhapcand; c++)
	{
		vector<int> final_mec(M,0);
		vector<int> allele(M,0);
		for (int i=0; i<M; i++)
		{
			vector<int> mec(K,0);
			for (int k=0; k<K; k++)
			{
				for (int j=0; j<=s; j++)
				{
					if (superread_mat[i][j]!=0)
					{
						if (superread_mat[i][j]!=hap[c][k][j])
							mec[k] = mec[k] + count_mat[i][j];
					}
				}
			}
			int val, ind;
			minIndVal(mec, val, ind);
			final_mec[i] = val;
			for (int j=0; j<=s; j++)
			{
				if (superread_mat[i][j]!=0)
					allele[i] = allele[i] + count_mat[i][j];
			}
		}
		MEC[c] = accumulate(final_mec.begin(),final_mec.end(),0);
		MECrate[c] = MEC[c]/double(accumulate(allele.begin(),allele.end(),0));
	}
}

// prune unlikely strains
void pruneHapSeq(int& numcand, vector<int>& remain, vector<vector<int> >& cand, int final_numhapcand, int K,vector<vector<vector<int>> >final_hap, int s)
{
	vector<int> candcheckind(numcand,0);
	for (int c=0; c<numcand; c++)
	{
		for (int i=0; i<final_numhapcand; i++)
		{
			if (candcheckind[c]==0)
			{
				for (int k=0; k<K; k++)
				{
					vector<int> mismatch;
					for (int j=0; j<=s; j++)
					{
						if (cand[c][j]!=final_hap[i][k][j])
							mismatch.push_back(c);
					}
					if (mismatch.size()==0)
					{
						candcheckind[c] = 1;
						break;
					}
				}
			}
			else
				break;
		}
	}
	remain.clear();
	int temp_numcand = numcand;
	numcand = 0;
	vector<vector<int> > temp_cand = cand;
	cand.clear();
	for (int c=0; c<temp_numcand; c++)
	{
		if (candcheckind[c]==1)
		{
			numcand++;
			remain.push_back(c);
			vector<int> temp_candvec(s+1);
			for (int j=0; j<=s; j++)
				temp_candvec[j] = temp_cand[c][j];
			cand.push_back(temp_candvec);
		}
	}
}

// candidates which can be potentially pruned further
void getPotentialPruneCand(vector<int>& potential_delete_order, vector<vector<int> > remain_comb, vector<double> final_pr, int final_numhapcand, int K)
{
	potential_delete_order.clear();
	for (int d=0; d<final_numhapcand-1; d++)
	{
		double val;
		int ind;
		minIndValdouble(final_pr, val, ind);
		vector<int> select(K);
		for (int k=0; k<K; k++)
			select[k] = remain_comb[ind][k];
		final_pr.erase(final_pr.begin()+ind);
		remain_comb.erase(remain_comb.begin()+ind);
		vector<int> del_flag(K,0);
		for (int i=0; i<K; i++)
		{
			int flag = 0;
			for (int j=0; j<remain_comb.size(); j++)
			{
				if (flag == 0)
				{
					for (int k=0; k<K; k++)
					{
						if (select[i]==remain_comb[j][k])
						{
							flag = 1;
							break;
						}
					}	
				}
			}
			if (flag == 1)
				del_flag[i] = 1;
		}
		for (int k=0; k<K; k++)
		{
			if (del_flag[k]==0)
				potential_delete_order.push_back(select[k]);
		}
	}
}

void pruneHapSeq2(int& numcand, vector<int>& remain, vector<vector<int> >& cand, int pre_numhapcand, int K,vector<vector<vector<int>> >final_hap, vector<vector<vector<int>> >pre_final_hap, int s, vector<double> final_pr)
{
	vector<int> candcheckind(numcand,0);
	for (int c=0; c<numcand; c++)
	{
		for (int i=0; i<pre_numhapcand; i++)
		{
			if (candcheckind[c]==0)
			{
				for (int k=0; k<K; k++)
				{
					vector<int> mismatch;
					for (int j=0; j<=s; j++)
					{
						if (cand[c][j]!=pre_final_hap[i][k][j])
							mismatch.push_back(c);
					}
					if (mismatch.size()==0)
					{
						candcheckind[c] = 1;
						break;
					}
				}
			}
			else
				break;
		}
		double val;
		int ind;
		maxIndValdouble(final_pr, val, ind);
		if (candcheckind[c]==0)
		{
			for (int k=0; k<K+1; k++)
			{
				vector<int> mismatch;
				for (int j=0; j<=s; j++)
				{
					if (cand[c][j]!=final_hap[ind][k][j])
						mismatch.push_back(c);
				}
				if (mismatch.size()==0)
				{
					candcheckind[c] = 1;
					break;
				}
			}

		}
	}
	remain.clear();
	int temp_numcand = numcand;
	numcand = 0;
	vector<vector<int> > temp_cand = cand;
	cand.clear();
	for (int c=0; c<temp_numcand; c++)
	{
		if (candcheckind[c]==1)
		{
			numcand++;
			remain.push_back(c);
			vector<int> temp_candvec(s+1);
			for (int j=0; j<=s; j++)
				temp_candvec[j] = temp_cand[c][j];
			cand.push_back(temp_candvec);
		}
	}
}

// sequential Baysian inference of quasispecies
void ML_QSR(vector<vector<int> > SNV_matrix, vector<vector<int> > SNV_freq, int nFrag, int nSNV, vector<vector<int> > superread_matrix, vector<vector<int> > count_matrix, vector<vector<double> > score_matrix, vector<int> nMember, double seq_err, vector<vector<int> >& haplo_snv, double eta1)
{
	double soft_eta1 = 0.5*eta1;
	int deletethres = round(nFrag*0.001);
	vector<int> keep;

	vector<vector<int> > superread_mat;
	vector<vector<int> > count_mat;
	vector<vector<double> > score_mat;
 	vector<int> nMem;

	vector<int> superread_vec(nSNV,0);
	vector<int> count_vec(nSNV,0);
	vector<double> score_vec(nSNV,0);

	for (int i=0; i<nMember.size(); i++)
	{
		if (nMember[i]>=deletethres)
		{
			keep.push_back(i);
			for (int j=0; j<nSNV; j++)
			{
				superread_vec[j] = superread_matrix[i][j];
				count_vec[j] = count_matrix[i][j];
				score_vec[j] = score_matrix[i][j];
			}
			superread_mat.push_back(superread_vec);
			count_mat.push_back(count_vec);
			score_mat.push_back(score_vec);
			nMem.push_back(nMember[i]);
		}
	}
	int M = nMem.size();
	int N = nSNV;
	vector<vector<double> > errtable(M,vector<double>(N));
	double logerr = log(1-seq_err);
	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (score_mat[i][j] < logerr)
				errtable[i][j] = exp(score_mat[i][j]);
			else
				errtable[i][j] = 1-seq_err;
		}
	}
	vector<vector<int> > pre_cand; // numcand by s
	vector<int> temp_precand;
	for (int i=0; i<4; i++)
	{
		vector<int> temp_precand(2,0);
		if (SNV_freq[0][i]>0)
		{
			temp_precand[0] = i+1;
			pre_cand.push_back(temp_precand);
		}
	}
	int pre_numcand = pre_cand.size();
	int K = pre_numcand;
	vector<int> remain(pre_numcand,0);
	for (int i=1; i<pre_numcand; i++)
		remain[i] = i;
	vector<int> potential_delete_order;
//	vector<int> track_K(N,0);               // to track MECvrate, MECimprate
//	vector<double> track_MECrate(N,0);      
//	vector<double> track_MECimprate(N,0);   

	for (int s=1; s<N; s++)
	{
       // cout << "At position "<< s << endl ;
		vector<int> permu;
		for (int i=0; i<4; i++)
		{
			if (SNV_freq[s][i]>0)
				permu.push_back(i+1);
		}
		vector<vector<int> > cand(pre_numcand*permu.size(),vector<int>(s+1));
		vector<int> preseq_index(pre_numcand*permu.size(),0);
		int numcand = 0;
		vector<int> vec_relNcov;
		vector<double> vec_relPR;
		getHapSeq(cand, numcand, M, superread_mat, pre_cand, pre_numcand, s, count_mat, permu, errtable, K, preseq_index, remain, vec_relNcov, vec_relPR);

		int numcomb_thres = maxNforKcomb(MAX_NUM_BRANCHES,K);
		while (numcand > numcomb_thres)    
		{ 
			if (potential_delete_order.size()>0)
			{   
				int idx = potential_delete_order[0];
				potential_delete_order.erase(potential_delete_order.begin());							
				vector<int> delseq;
				for (int i=0; i<preseq_index.size(); i++)
				{ 
					if (preseq_index[i]==idx)
						delseq.push_back(i);
				}
				numcand = numcand-delseq.size();
				
				for (int i=delseq.size()-1; i>=0; i=i-1 )
				{
					cand.erase(cand.begin()+delseq[i]);
					preseq_index.erase(preseq_index.begin()+delseq[i]);
					vec_relNcov.erase(vec_relNcov.begin()+delseq[i]);
					vec_relPR.erase(vec_relPR.begin()+delseq[i]);
				}				

			}	
			else
			{   
				int minval, minind;
				minIndVal(vec_relNcov, minval, minind);
				vector<double> relPR_list;
				vector<int> indlist;
				for (int i=0; i<vec_relNcov.size(); i++)
				{
					if (vec_relNcov[i]==minval)
					{
						indlist.push_back(i);
						relPR_list.push_back(vec_relPR[i]);
					}
				}
				double val;
				int ind;
				minIndValdouble(relPR_list, val, ind);
				int del_ind = indlist[ind];		
				cand.erase(cand.begin()+del_ind);
				vec_relNcov.erase(vec_relNcov.begin()+del_ind);
				vec_relPR.erase(vec_relPR.begin()+del_ind);
				numcand = numcand-1;			
			}
		}
		vector<double> mecrate;
		vector<double> mecimprate;
		vector<int> mecK;
		
		int final_numhapcand = 0, pre_numhapcand;
		vector<double> final_pr, pre_final_pr;
		vector<vector<vector<int> > > final_hap, pre_final_hap;
		vector<vector<int> > remain_comb, pre_remain_comb;
		for (int Kiter = 0; Kiter<20; Kiter++)  // Kiter can be increased
		{
			vector<vector<int> > comb;
			if (numcand == K)
			{
				vector<int> combvec;	
				for (int i=0; i<K; i++)
					combvec.push_back(i);	
				comb.push_back(combvec);
			}
			else
				combnk(numcand,K,comb);
			int numhapcand = comb.size();
			vector<vector<vector<int> > > hap; // numhapcand by K by s+1
			for (int i=0; i<numhapcand; i++)
			{
				vector<int> sel;
				vector<vector<int> > hap_mat(K,vector<int>(s+1));
				for (int j=0; j<K; j++)
				{
					sel.push_back(comb[i][j]);
					for (int k=0; k<s+1; k++)
					{	
						hap_mat[j][k] = cand[comb[i][j]][k];
					}
				}
				hap.push_back(hap_mat);
		
				
			}
			if (numhapcand == 1) ///////////////////////check values,index.. 
			{
				vector<vector<int> > hap_mat(K,vector<int>(s+1));
				hap.push_back(hap_mat);
			} // hap : 2 by K by s+1.. hap[2][][] is filled with 0 for 3 dim

			vector<double> logPR(numhapcand,0);
			getLogLikeli(logPR, numhapcand, M, K, s, superread_mat, hap, errtable);
			double maxPR;
			int maxPRind;
			maxIndValdouble(logPR, maxPR, maxPRind);

			vector<int> MEC(numhapcand,0);
			vector<double> MECrate(numhapcand,0);
			getMECrate(MEC, MECrate, numhapcand, M, K, s, superread_mat, hap, count_mat);

			mecK.push_back(MEC[maxPRind]);
			mecrate.push_back(MECrate[maxPRind]);
			
			final_numhapcand = 0;
			final_pr.clear();
			final_hap.clear(); 
			remain_comb.clear();
			for (int c=0; c<numhapcand; c++)
			{
				if ( exp(logPR[c]-maxPR) > delta1)
				{
					final_numhapcand++;

					final_pr.push_back(logPR[c]);
					vector<vector<int> > temp_hap(K,vector<int>(s+1));
					for (int k=0;k<K;k++)
					{
						for (int j=0; j<=s; j++)
							temp_hap[k][j] = hap[c][k][j];
					}
					final_hap.push_back(temp_hap);
					vector<int> temp_comb(K);
					for (int k=0; k<K; k++)
						temp_comb[k] = comb[c][k];
					remain_comb.push_back(temp_comb);
				}
			}

			
			if (Kiter == 0)
			{
				if (mecrate[Kiter] < seq_err*2 || K==numcand)
				{
					pruneHapSeq(numcand, remain, cand, final_numhapcand, K, final_hap, s);
					
					getPotentialPruneCand(potential_delete_order, remain_comb, final_pr, final_numhapcand, K);
				//	track_K[s] = K;
				//	track_MECrate[s] = mecrate[Kiter];
					break;
				}
				else
				{
					K++;	
					pre_final_pr = final_pr;
					pre_final_hap = final_hap;
					pre_numhapcand = final_numhapcand;
					pre_remain_comb = remain_comb;
					continue;
				}
			}

			if (Kiter > 0)
			{
				double MECimprate = (mecK[Kiter-1]-mecK[Kiter])/double(mecK[Kiter-1]);
				mecimprate.push_back(MECimprate);
						if (MECimprate >= eta1)
				{
					if (MECimprate == 1 || mecrate[Kiter] < seq_err*2 || K == numcand)
					{
						pruneHapSeq(numcand, remain, cand, final_numhapcand, K, final_hap, s);
					
						getPotentialPruneCand(potential_delete_order, remain_comb, final_pr, final_numhapcand, K);
					//	track_K[s] = K;
					//	track_MECrate[s] = mecrate[Kiter];
					//	track_MECimprate[s] = MECimprate;
						break;
					}	
				}
				else
				{
					if (K-2<=1)
					{
						K = K-1;
						pruneHapSeq(numcand, remain, cand, pre_numhapcand, K, pre_final_hap, s);

						getPotentialPruneCand(potential_delete_order, pre_remain_comb, pre_final_pr, pre_numhapcand, K);
						final_pr = pre_final_pr;
						final_hap = pre_final_hap;
	
					//	track_K[s] = K;
					//	track_MECrate[s] = mecrate[Kiter-1];
						break;			
					}

					if (Kiter < 2)
					{ 
						K = K-1;
						int pre_K = K;
						while (MECimprate < eta1)
						{
							K = K-1;
							comb.clear();
							if (numcand == K)
							{
								vector<int> combvec;	
								for (int i=0; i<K; i++)
									combvec.push_back(i);	
								comb.push_back(combvec);
							}
							else
								combnk(numcand,K,comb);
							numhapcand = comb.size();
		        				hap.clear(); // numhapcand by K by s+1
							for (int i=0; i<numhapcand; i++)
							{
								vector<int> sel;
								vector<vector<int> > hap_mat(K,vector<int>(s+1));
								for (int j=0; j<K; j++)
								{
									sel.push_back(comb[i][j]);
									for (int k=0; k<s+1; k++)
									{	
										hap_mat[j][k] = cand[comb[i][j]][k];
									}
								}
								hap.push_back(hap_mat);
							}
							if (numhapcand == 1) 
							{
								vector<vector<int> > hap_mat(K,vector<int>(s+1));
								hap.push_back(hap_mat);
							} // hap : 2 by K by s+1.. hap[2][][] is filled with 0 for 3 dim

							logPR.clear();
							logPR.assign(numhapcand,0);
							getLogLikeli(logPR, numhapcand, M, K, s, superread_mat, hap, errtable);
							maxIndValdouble(logPR, maxPR, maxPRind);
							MEC.clear();
							MEC.assign(numhapcand,0);
							MECrate.clear();
							MECrate.assign(numhapcand,0);
							getMECrate(MEC, MECrate, numhapcand, M, K, s, superread_mat, hap, count_mat);
							mecK.insert (mecK.begin(),MEC[maxPRind]);
							mecrate.insert (mecrate.begin(),MECrate[maxPRind]);
							MECimprate = (mecK[0]-mecK[1])/double(mecK[0]);
							mecimprate.insert (mecimprate.begin(),MECimprate);		
							final_numhapcand = 0;
							final_pr.clear();
							final_hap.clear(); 
							remain_comb.clear();
							for (int c=0; c<numhapcand; c++)
							{
								if ( exp(logPR[c]-maxPR) > delta1)
								{
									final_numhapcand++;

									final_pr.push_back(logPR[c]);
									vector<vector<int> > temp_hap(K,vector<int>(s+1));
									for (int k=0;k<K;k++)
									{
										for (int j=0; j<=s; j++)
											temp_hap[k][j] = hap[c][k][j];
									}
									final_hap.push_back(temp_hap);
									vector<int> temp_comb(K);
									for (int k=0; k<K; k++)
										temp_comb[k] = comb[c][k];
									remain_comb.push_back(temp_comb);
								}
							}
							if (MECimprate < eta1)
							{
								pre_K = K;
								pre_final_pr = final_pr;
								pre_final_hap = final_hap;
								pre_numhapcand = final_numhapcand;
								pre_remain_comb = remain_comb;
							}
						}
						K = pre_K;

						if (mecimprate[1] > soft_eta1)
						{
							K++;
							comb.clear();
							if (numcand == K)
							{
								vector<int> combvec;	
								for (int i=0; i<K; i++)
									combvec.push_back(i);	
								comb.push_back(combvec);
							}
							else
								combnk(numcand,K,comb);
							numhapcand = comb.size();
		        				hap.clear(); // numhapcand by K by s+1
							for (int i=0; i<numhapcand; i++)
							{
								vector<int> sel;
								vector<vector<int> > hap_mat(K,vector<int>(s+1));
								for (int j=0; j<K; j++)
								{
									sel.push_back(comb[i][j]);
									for (int k=0; k<s+1; k++)
									{	
										hap_mat[j][k] = cand[comb[i][j]][k];
									}
								}
								hap.push_back(hap_mat);
							}
							if (numhapcand == 1)  
							{
								vector<vector<int> > hap_mat(K,vector<int>(s+1));
								hap.push_back(hap_mat);
							} // hap : 2 by K by s+1.. hap[2][][] is filled with 0 for 3 dim

							logPR.clear();
							logPR.assign(numhapcand,0);
							getLogLikeli(logPR, numhapcand, M, K, s, superread_mat, hap, errtable);
							maxIndValdouble(logPR, maxPR, maxPRind);
							final_numhapcand = 0;
							final_pr.clear();
							final_hap.clear(); 
							remain_comb.clear();
							for (int c=0; c<numhapcand; c++)
							{
								if ( exp(logPR[c]-maxPR) > delta1)
								{
									final_numhapcand++;

									final_pr.push_back(logPR[c]);
									vector<vector<int> > temp_hap(K,vector<int>(s+1));
									for (int k=0;k<K;k++)
									{
										for (int j=0; j<=s; j++)
											temp_hap[k][j] = hap[c][k][j];
									}
									final_hap.push_back(temp_hap);
									vector<int> temp_comb(K);
									for (int k=0; k<K; k++)
										temp_comb[k] = comb[c][k];
									remain_comb.push_back(temp_comb);
								}
							}

							K = pre_K;

							pruneHapSeq2(numcand, remain, cand, pre_numhapcand, K, final_hap, pre_final_hap, s, final_pr);

							final_pr = pre_final_pr;
						}
						else
						{
							K = pre_K;
							final_pr = pre_final_pr;

							pruneHapSeq(numcand, remain, cand, pre_numhapcand, K, pre_final_hap, s);

						}
						getPotentialPruneCand(potential_delete_order, pre_remain_comb, pre_final_pr, pre_numhapcand, K);

						final_pr = pre_final_pr;
						final_hap = pre_final_hap;
	
					//	track_K[s] = K;
					//	track_MECrate[s] = mecrate[1];
					//	track_MECimprate[s] = MECimprate;
						break;	
					}
					else 
					{	
						K = K-1;
						pruneHapSeq(numcand, remain, cand, pre_numhapcand, K, pre_final_hap, s);

						getPotentialPruneCand(potential_delete_order, pre_remain_comb, pre_final_pr, pre_numhapcand, K);
						final_pr = pre_final_pr;
						final_hap = pre_final_hap;
	
					//	track_K[s] = K;
					//	track_MECrate[s] = mecrate[Kiter-1];
					//	track_MECimprate[s] = mecimprate[Kiter-2];
						break;	
					}

				}
			}
			K++;
			pre_final_pr = final_pr;
			pre_final_hap = final_hap;

			pre_numhapcand = final_numhapcand;
			pre_remain_comb = remain_comb;						
		}
		pre_cand = cand;
		pre_numcand = numcand;

		if (s == N-1)
		{
			double val;
			int ind;
			maxIndValdouble(final_pr, val, ind);
			for (int k=0; k<K; k++)
			{
				vector<int> haplo_seq(N);
				for (int s=0; s<N; s++)
					haplo_seq[s] = final_hap[ind][k][s];
				haplo_snv.push_back(haplo_seq);
			}
		}
	}

}

// estimate relative frequencies of strains
void estimateViralFreq(vector<double>& viral_freq, vector<vector<int> > viral_quasi, vector<vector<int>> haplo_snv, vector<vector<int>> Read_matrix, vector<vector<int> > superread_matrix, vector<int> nMember, int nGroup, vector<vector<int> > group_matrix, int num_strain, int nSNV, int gene_length, int nReads)
{
	vector<int> viral_abun(num_strain,0);
	vector<int> assigned_ind(nReads,0);
	for (int i=0; i<nGroup; i++)
	{
		vector<int> dist(num_strain,0);
		for (int k=0; k<num_strain; k++)
		{
			int HD =0;
			for (int j=0; j<nSNV; j++)
			{	
				if (superread_matrix[i][j]!=haplo_snv[k][j])	
					HD++; 
			}
			dist[k] = HD;
		}
		int val, ind;
		minIndVal(dist, val, ind);
		viral_abun[ind] = viral_abun[ind] + nMember[i];

		for (int j=0; j<nMember[i]; j++)
		{
			assigned_ind[group_matrix[i][j]] = 1;
		}
	}
	vector<int> unassigned_list;
	for (int i=0; i<nReads; i++)
	{	
		if (assigned_ind[i]==0)
			unassigned_list.push_back(i);
	}
	
	for (int i=0; i<unassigned_list.size(); i++)
	{
		vector<int> dist(num_strain,0);
		for (int k=0; k<num_strain; k++)
		{
			int HD =0;
			for (int j=0; j<gene_length; j++)
			{	
				if (Read_matrix[unassigned_list[i]][j]!=viral_quasi[k][j])
					HD++; 
			}
			dist[k] = HD;
		}
		int val, ind;
		minIndVal(dist, val, ind);
		viral_abun[ind] = viral_abun[ind] + 1;
	}
	for (int i=0; i<num_strain; i++)
		viral_freq[i] = viral_abun[i]/double(nReads);
}

void sortFreq(vector<double> viral_freq, int num_strain, vector<int>& index)
{
	int j, P;
	int index_temp;
	double Tmp;
	for (P = 1; P < num_strain; P++)
	{
		Tmp = viral_freq[P];
		index_temp = P;
		for (j = P; j > 0 && viral_freq[j-1] < Tmp; j--)
		{
			viral_freq[j] = viral_freq[j-1];
			index[j] = index[j-1];
		}
		viral_freq[j] = Tmp;
		index[j] = index_temp;
	}
}

void makeoutput(vector<vector<int> > viral_quasi, vector<double> viral_freq, int num_strain, int gene_length, string zonename)
{
	vector<int> index(num_strain,0);
	sortFreq(viral_freq, num_strain, index);
	/*for (int i=0; i<num_strain; i++)
		cout << index[i] << " ";
	cout << endl;
     */
    
    cout << "Estimated population size: " << num_strain << endl;
   // cout << "Estimated frequencies "<< endl ;
    for (int i=0; i<num_strain; i++)
    {
        cout << "Frequency of strain " << i+1 << ": " << viral_freq[index[i]]<< endl;
    }
    cout << endl;

	std::ofstream writefile1;
	//std::ofstream writefile2;
	//std::ofstream writefile3;

	std::string name = zonename;

	writefile1.open(name+"_ViralSeq.txt");
//	writefile2.open(n ame+"_Seq.txt");
//	writefile3.open(name+"_Freq.txt");


	for (int i=0; i<num_strain; i++)
	{
		writefile1 << "Viral Quasispecies - strain"<< i+1 << "_freq: " << viral_freq[index[i]];
		writefile1 << "\n";
		for (int j=0; j<gene_length; j++)
		{
			switch (viral_quasi[index[i]][j])
			{	
				case 1:
					writefile1 << "A";
					break;
				case 2:
					writefile1 << "C";
					break;
				case 3:
					writefile1 << "G";
					break;
				case 4:
					writefile1 << "T";
					break;
				default:
					writefile1 << "N";
					break;
			}
		}
		writefile1 << "\n";
	}

/*	for (int i=0; i<num_strain; i++)
	{
		for (int j=0; j<gene_length; j++)
		{
			switch (viral_quasi[index[i]][j])
			{	
				case 1:
					writefile2 << "A";
					break;
				case 2:
					writefile2 << "C";
					break;
				case 3:
					writefile2 << "G";
					break;
				case 4:
					writefile2 << "T";
					break;
				default:
					writefile2 << "N";
					break;
			}
		}
		writefile2 << "\n";
	}

	for (int i=0; i<num_strain; i++)
		writefile3 << viral_freq[index[i]] << " ";
*/
	writefile1.close();
//	writefile2.close();
//	writefile3.close();		
}

int main(int argc, char* argv[]) {
	if(argc < 2)
	{
    		cout <<"usage: aBayesQR <config.txt>" << endl;
    		return 0;
  	}
 
  	srand(time(NULL));
 
  	string line, line_stats, line_ID;
  	string tok = ":";
  	vector<string> tokens, tokens2, tokens_as;
  
  	string confStr;
 
	if(argc == 2) confStr = argv[1];

	ifstream infConf(confStr.c_str(), ios::in);
  	vector<string> arg_buffer;
	int count =0;
	string arg_buffer_sub;
	
	cout << endl;
	while (!infConf.eof() )
	{
		getline(infConf,line,'\n');
	    	if( line.length() > 0 ) 
		{
	      		size_t pos = line.find(":");
	      		cout << count++ << ". "<< line << endl;
	      		arg_buffer_sub = line.substr(pos+2);
	      		arg_buffer.push_back(arg_buffer_sub);
	      		arg_buffer_sub.clear();
	    	}
  	}
  	infConf.close();
	
	string cons;
	vector<string> FASTAreads;
	int pairedend;
	double SNV_thres;
	double overlap_ratio_start = 0.9, overlap_ratio_end = 0.1;
	int  reconstruction_start,  reconstruction_end;
	double  min_qual;
	int min_length, max_insertln;
	double  max_gap_fraction = 0.05, min_align_score_fraction = 0.35;
	double seq_err;
	double eta1;
	
	string zonename;
	
	if(count < 3){
	    cout <<"problem with config file...please check" << endl;
	    return 0;
	}
	cons =  arg_buffer[0];
	FASTAreads.push_back(arg_buffer[1]);
	pairedend = atoi(arg_buffer[2].c_str());
	SNV_thres = atof(arg_buffer[3].c_str());

	if(count > 4){
	 //   overlap_ratio_start  = atof(arg_buffer[3].c_str());
	  //  if(count > 4){
	  //    overlap_ratio_end = atof(arg_buffer[4].c_str());
	  //    if(count > 5){
		reconstruction_start = atoi(arg_buffer[4].c_str());
		if(count > 5){
		  reconstruction_end = atoi(arg_buffer[5].c_str());
		  if(count > 6){
		    min_qual = atoi(arg_buffer[6].c_str());
		    if(count > 7){
		      min_length =  atoi(arg_buffer[7].c_str());
		      if(count > 8){
		      	max_insertln = atoi(arg_buffer[8].c_str());
			//if(count > 10){
			//  max_gap_fraction = atof(arg_buffer[10].c_str());
			//  if(count > 11){
			//  min_align_score_fraction = atof(arg_buffer[11].c_str());
				if(count >9){
				zonename = arg_buffer[9];
					if (count > 10){
					seq_err = atof(arg_buffer[10].c_str());
						if (count > 11){
							eta1 = atof(arg_buffer[11].c_str());
						}
					}
				}
			  
			//}
		//	}
		      }
		    }
		  }
		//}
	    //}
	  }
	}
	seq_err = seq_err*0.01;	
	
	clock_t start, end;
	double cpu_time_used;
	start = clock();

	vector<vector<int> > Read_matrix;	  
	int total_count = 0;
	char gap_quality = '*'; //'I'; 	 
	double mean_length = 0.;
	int gene_length  = reconstruction_end-reconstruction_start+1;
	
	int error_flag = 0;
	if (pairedend == 1)
  		error_flag = parseSAMpaired(FASTAreads[0],  max_gap_fraction, min_align_score_fraction, min_qual,  min_length, max_insertln, gap_quality,  mean_length, Read_matrix, reconstruction_start, reconstruction_end, total_count, gene_length);
	else
		error_flag = parseSAM(FASTAreads[0],  max_gap_fraction, min_align_score_fraction, min_qual,  min_length, max_insertln, gap_quality,  mean_length, Read_matrix, reconstruction_start, reconstruction_end, total_count, gene_length);
	
 	if(error_flag>0)
    		return 1;
	int nReads = Read_matrix.size();

	for (int i = 0; i< Read_matrix.size(); i++)
	{
	        if (Read_matrix[i].size() != gene_length)
			cout << "error!!!"<<endl;
	}
	if (pairedend == 1)
		cout <<  endl << "After parsing " << total_count << " sequencing reads in file " <<FASTAreads[0]<< ", there are "<<nReads<< " paired_end reads(mean lengths "<< mean_length/Read_matrix.size() << ") covering regions "<< reconstruction_start << "-" << reconstruction_end <<"."<< endl;
	else
		cout <<  endl << "After parsing " << total_count << " sequencing reads in file " <<FASTAreads[0]<< ", there are "<<nReads<< " reads(mean lengths "<< mean_length/Read_matrix.size() << ") covering regions "<< reconstruction_start << "-" << reconstruction_end <<"."<< endl;
	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for SAM parsing: "  << cpu_time_used << endl<<endl;

	start = clock();		
	vector<vector<int> > Allele_freq; // nReads by 5
	vector<vector<int> > SNV_matrix; // nFrag by nSNV
	vector<vector<int> > SNV_freq; // nSNV by 5
	vector<int> Homo_seq;
	vector<int> SNV_pos;
	int nSNV;
	vector<int> deleted_reads_list;	
	callSNV(Read_matrix, Allele_freq, SNV_matrix, SNV_freq, Homo_seq, SNV_pos, nReads, gene_length, SNV_thres, nSNV,deleted_reads_list);
	int nFrag = SNV_matrix.size();	
	cout << "After calling SNVs from " << gene_length << " bases in regions between " << reconstruction_start << " and " << reconstruction_end << ", " << nSNV<< " SNVs are detected." << endl;
	cout << "After correcting error, "<< nFrag <<" fragments are used for quasi-species reconstruction." << endl;
	cout << "reduced number of fragment:" << deleted_reads_list.size() << endl;
	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for SNV calling: "  << cpu_time_used << endl << endl;
	start = clock();
	vector<vector<int> > group_matrix;
    	vector<vector<int> > superread_matrix;
    	vector<vector<int> > count_matrix;
    	vector<vector<double> > score_matrix;
    	vector<int> nMember;
	int nGroup = 0;
	aggromerative_clustering(SNV_matrix, SNV_freq, nFrag, nSNV, overlap_ratio_start, overlap_ratio_end, zonename, min_length, group_matrix, superread_matrix, count_matrix, score_matrix, nMember, nGroup);
	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for agglomerative clustering for zone " << zonename <<" : "  << cpu_time_used << endl << endl;
    	cout << "total number of group: " << nMember.size() << endl;
    
    	start = clock();
    cout << endl<< "Quasispecies reconstruction starts with " << nMember.size()<< " superreads covering "<< nSNV<< " SNVs"  << endl;
	vector<vector<int>> haplo_snv;
    	ML_QSR(SNV_matrix, SNV_freq, nFrag, nSNV, superread_matrix, count_matrix, score_matrix, nMember, seq_err, haplo_snv, eta1);
	end = clock();  
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        cout << endl << "CPU time for sequential Bayesian inference for zone " << zonename <<" : "  << cpu_time_used << endl << endl;

	int num_strain = haplo_snv.size();
	vector<vector<int> > viral_quasi;
	for (int i=0; i<num_strain; i++)
		viral_quasi.push_back(Homo_seq);
	for (int i=0; i<num_strain; i++)
	{
		for (int j=0; j<nSNV; j++)
		{	
			viral_quasi[i][SNV_pos[j]] = haplo_snv[i][j];
		}
	}
	/*cout << "Viral Quasispecies: "   << endl;
	for (int i=0; i<num_strain; i++)
	{
		for (int j=0; j<gene_length; j++)
		{	
			cout << viral_quasi[i][j];
		}
		cout << endl;
	}*/

	vector<double> viral_freq(num_strain,0);
	estimateViralFreq(viral_freq, viral_quasi, haplo_snv, Read_matrix, superread_matrix, nMember, nGroup, group_matrix, num_strain, nSNV, gene_length, nReads);
	
	/*cout << "Viral frequency : "   << endl;
	for (int i=0; i<num_strain; i++)
		cout << viral_freq[i] << " ";
	cout << endl;*/

	makeoutput(viral_quasi, viral_freq, num_strain, gene_length, zonename);
}
