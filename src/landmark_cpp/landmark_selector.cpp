#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<string>
#include<sstream>
#include<cstdlib>
#include<omp.h>
#include<cmath>
#include<cstdio>
#include<cfloat>
#include "random.h"

using namespace std;

struct IndexValue{
	int index;
	float value;
	IndexValue(int ind, float val){
		index = ind;
		value = val;
	}
};

bool sort_ascending(IndexValue i,IndexValue j) { return (i.value < j.value); }
// 
// class InputParser{
//     public:
//         InputParser (int &argc, char **argv){
//             for (int i=1; i < argc; ++i)
//                 this->tokens.push_back(string(argv[i]));
//         }
// 
//         const string& getCmdOption(const string &option) const{
//             vector<string>::const_iterator itr;
//             itr =  find(this->tokens.begin(), this->tokens.end(), option);
//             if (itr != this->tokens.end() && ++itr != this->tokens.end()){
//                 return *itr;
//             }
//             static const string empty_string("");
//             return empty_string;
//         }
// 
//         bool cmdOptionExists(const string &option) const{
//             return find(this->tokens.begin(), this->tokens.end(), option)
//                    != this->tokens.end();
//         }
//     private:
//         vector <string> tokens;
// };
// 
// 
// void parse_string_2_int(string s, vector<int> &res)
// {
// 	istringstream str_stream(s);
// 	while(str_stream)
// 	{
// 		string tmp;
// 		if (!getline(str_stream, tmp, ',')) break;
// 		res.push_back(stoi(tmp));
// 	}
// }
// 
// void parse_string_2_float(string s, vector<float> &res)
// {
// 	istringstream str_stream(s);
// 	while(str_stream)
// 	{
// 		string tmp;
// 		if (!getline(str_stream, tmp, ',')) break;
// 		res.push_back(stof(tmp));
// 	}
// }

// [start, end)
float compute_distance(const vector<float> &data, const vector<int> &indices,
		int start_i, int end_i, int start_j, int end_j)
{
	float dist = 0.0;
	int pos_i = start_i;
	int pos_j = start_j;
	while (pos_i < end_i && pos_j < end_j){
		if (indices[pos_i] == indices[pos_j])
		{
			dist += (data[pos_i] - data[pos_j]) * (data[pos_i] - data[pos_j]);
			++pos_i;
			++pos_j;
		}
		else if(indices[pos_i] < indices[pos_j])
		{
			dist += data[pos_i] * data[pos_i];
			++pos_i;
		}else
		{
			dist += data[pos_j] * data[pos_j];
			++pos_j;
		}
	}

	if (pos_i < end_i)
	{
		for(int t=pos_i; t<end_i;t++) dist += data[t] * data[t];
	}else{
		for(int t=pos_j; t<end_j;t++) dist += data[t] * data[t];
	}

	return sqrt(dist);
}

void select_landmarks(const vector<float> &data, const vector<int> &indices,
		const vector<int> &indptr, int n, int dim, int count,
		vector<float> &dist, vector<int> &assign, vector<int> &flag){

	// select the first landmark
	int landmark = random() % n;

#pragma omp parallel for
	for(int i=0;i<n;i++){
//		int tid = omp_get_thread_num();
//		printf("example %d,thread %d\n",i, tid);

		if (flag[i]==1) continue;
		if (i==landmark)
		{
			dist[i] = 0.0;
			assign[i] = landmark;
			continue;
		}
		float dist_i_landmark = compute_distance(data, indices, indptr[i], indptr[i+1], indptr[landmark], indptr[landmark+1]);
		if (dist_i_landmark < dist[i]){
			dist[i] = dist_i_landmark;
			assign[i] = landmark;
		}
	}
	flag[landmark] = 1;
	// cout<<landmark<<endl;

	// select the rest of landmarks
	for(int t=1;t<count;t++){

		// sort distance in ascending order for non-landmarks
		vector<IndexValue> indval;
		for(int i=0;i<n;i++){
			if (flag[i]==1) continue;
			IndexValue tmp(i,dist[i]);
			indval.push_back(tmp);
		}
		sort(indval.begin(), indval.end(), sort_ascending);

		int total = indval.size();
		int half = total/2;
		int number = random() % half + (total-half);
		landmark = indval[number].index;

//		for (size_t ii = 0;ii<indval.size();ii++)
//		{
//			cout<<indval[ii].index << ","<<indval[ii].value<<endl;
//		}
//		cout<<"landmark="<<landmark<<",number="<<number<<endl;
//		if (t > 2) exit(0);

		// update distances
#pragma omp parallel for
		for(int i=0;i<n;i++){
			if (flag[i]==1) continue;
			if (i==landmark)
			{
				dist[i] = 0.0;
				assign[i] = landmark;
				continue;
			}
			float dist_i_landmark = compute_distance(data, indices, indptr[i], indptr[i+1], indptr[landmark], indptr[landmark+1]);
			if (dist_i_landmark < dist[i]){
				dist[i] = dist_i_landmark;
				assign[i] = landmark;
			}
		}
		flag[landmark] = 1;
		// cout<<landmark<<endl;
	}
}
// 
// int test2()
// {
// 	int total = omp_get_max_threads();
// 	cout<<total;
// 
// 	int i;
// 	int numthreads = total;
// #pragma omp parallel for default(none) num_threads(numthreads) private(i)
// 	for (i = 0; i < 100; i++)
// 	{
// 		int tid = omp_get_thread_num();
// 		printf("Hello world from omp thread %d\n", tid);
// 	}
// 
// 	return 0;
// }
// 
// int main(int argc, char**argv){
// 
// 	InputParser parser(argc, argv);
// 	if(parser.cmdOptionExists("-h")){
// 		cout<<"command format:"<<endl;
// 		cout<<"landmark_selector -p indptr-file -d data-file -i indices-file -c cores -n landmarks -o output-dir"<<endl;
// 	}
// 
// 	string indptr_file = parser.getCmdOption("-p");
// 	if(indptr_file.empty()){
// 		cout<<"missing indptr file"<<endl;
// 		return 0;
// 	}
// 	string data_file = parser.getCmdOption("-d");
// 	if(data_file.empty()){
// 		cout<<"missing data file"<<endl;
// 		return 0;
// 	}
// 	string indices_file = parser.getCmdOption("-i");
// 	if(indices_file.empty()){
// 		cout<<"missing indices file"<<endl;
// 		return 0;
// 	}
// 
// 	string str_count = parser.getCmdOption("-n");
// 	if (str_count.empty()){
// 		return 0;
// 	}
// 	int count = stoi(str_count);
// 
// 	string str_num_proc = parser.getCmdOption("-c");
// 	if (str_num_proc.empty()){
// 		return 0;
// 	}
// 	int num_proc = stoi(str_num_proc);
// 
// 	string output_dir = parser.getCmdOption("-o");
// 	if (output_dir.empty()){
// 		cout<<"missing output directory"<<endl;
// 		return 0;
// 	}
// 
// //	int count = 10;
// //	int num_proc = 8;
// 	omp_set_dynamic(0);     // Explicitly disable dynamic teams
// 	omp_set_num_threads(num_proc);
// 
// 	double t1 = omp_get_wtime();
// 	string line;
// 
// 	cout<<"loading "<<indptr_file<<endl;
// 	ifstream findptr(indptr_file,ifstream::in);
// 	getline(findptr,line);
// 	vector<int> d_n;
// 	parse_string_2_int(line, d_n);
// 	int dim = d_n[0];
// 	int n = d_n[1];
// 	cout<<"n="<<n<<",dim="<<dim<<endl;
// 
// 	vector<int> indptr;
// 	getline(findptr,line);
// 	parse_string_2_int(line, indptr);
// 	findptr.close();
// 
// 	cout<<"loading "<<data_file<<endl;
// 	ifstream fdata(data_file,ifstream::in);
// 	getline(fdata,line);
// 	vector<float> data;
// 	parse_string_2_float(line, data);
// 	fdata.close();
// 
// 	cout<<"loading "<<indices_file<<endl;
// 	ifstream findices(indices_file,ifstream::in);
// 	getline(findices,line);
// 	vector<int> indices;
// 	parse_string_2_int(line, indices);
// 	findices.close();
// 	t1 = omp_get_wtime() - t1;
// 
// 	double t = omp_get_wtime();
// 	vector<float> dist(n, FLT_MAX); // will this cause a lot memory footprint? 
// 	vector<int> assign(n, -1);
// 	// 0 for non-landmark, 1 for landmark
// 	vector<int> flag(n, 0);
// 	select_landmarks(data, indices, indptr, n, dim, count, dist, assign, flag);
// 	t = omp_get_wtime() - t;
// 
// 	// output files
// 	double t2 = omp_get_wtime();
// 	string landmark_file = output_dir+"/landmark.txt";
// 	ofstream out_landmark(landmark_file, ofstream::out);
// 	for(int i = 0;i<n;i++){
// 		if (flag[i] == 1){
// 			out_landmark<<i;
// 			for(int pos=indptr[i];pos<indptr[i+1];pos++){
// 				out_landmark<<" "<<indices[pos]<<":"<<data[pos];
// 			}
// 			out_landmark<<"\n";
// 		}
// 	}
// 	string assign_dist = output_dir+"/assign_dist.txt";
// 	ofstream out_dist(assign_dist, ofstream::out);
// 	for(int i=0;i<n;i++){
// 		out_dist<<i<<" "<<assign[i]<<" "<<dist[i]<<"\n";
// 	}
// 	t2 = omp_get_wtime() - t2;
// 
// 	cout<<"loading data takes "<<t1<<" seconds"<<endl;
// 	cout<<"processing takes "<<t<<" seconds"<<endl;
// 	cout<<"output data takes "<<t2<<" seconds"<<endl;
// 	printf ("It took me %f seconds in total.\n",t1+t+t2);
// 
// 	return 0;
// }
