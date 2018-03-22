#include <Rcpp.h>
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
#include "./random.h"
#include<omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
  */


struct IndexValue{
  int index;
  float value;
  IndexValue(int ind, float val){
    index = ind;
    value = val;
  }
};

bool sort_ascending(IndexValue i,IndexValue j) { return (i.value < j.value); }

class InputParser{
public:
  InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i)
      this->tokens.push_back(string(argv[i]));
  }
  
  const string& getCmdOption(const string &option) const{
    vector<string>::const_iterator itr;
    itr =  find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
      return *itr;
    }
    static const string empty_string("");
    return empty_string;
  }
  
  bool cmdOptionExists(const string &option) const{
    return find(this->tokens.begin(), this->tokens.end(), option)
    != this->tokens.end();
  }
private:
  vector <string> tokens;
};


void parse_string_2_int(string s, vector<int> &res)
{
  istringstream str_stream(s);
  while(str_stream)
  {
    string tmp;
    if (!getline(str_stream, tmp, ',')) break;
    res.push_back(stoi(tmp));
  }
}

void parse_string_2_float(string s, vector<float> &res)
{
  istringstream str_stream(s);
  while(str_stream)
  {
    string tmp;
    if (!getline(str_stream, tmp, ',')) break;
    res.push_back(stof(tmp));
  }
}

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

void select_landmarks_cpp(const vector<float> &data, const vector<int> &indices,
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

int test2()
{
  int total = omp_get_max_threads();
  cout<<total;
  
  int i;
  int numthreads = total;
#pragma omp parallel for default(none) num_threads(numthreads) private(i)
  for (i = 0; i < 100; i++)
  {
    int tid = omp_get_thread_num();
    printf("Hello world from omp thread %d\n", tid);
  }
  
  return 0;
}

int main(int argc, char**argv){
  
  InputParser parser(argc, argv);
  if(parser.cmdOptionExists("-h")){
    cout<<"command format:"<<endl;
    cout<<"landmark_selector -p indptr-file -d data-file -i indices-file -c cores -n landmarks -o output-dir"<<endl;
  }
  
  string indptr_file = parser.getCmdOption("-p");
  if(indptr_file.empty()){
    cout<<"missing indptr file"<<endl;
    return 0;
  }
  string data_file = parser.getCmdOption("-d");
  if(data_file.empty()){
    cout<<"missing data file"<<endl;
    return 0;
  }
  string indices_file = parser.getCmdOption("-i");
  if(indices_file.empty()){
    cout<<"missing indices file"<<endl;
    return 0;
  }
  
  string str_count = parser.getCmdOption("-n");
  if (str_count.empty()){
    return 0;
  }
  int count = stoi(str_count);
  
  string str_num_proc = parser.getCmdOption("-c");
  if (str_num_proc.empty()){
    return 0;
  }
  int num_proc = stoi(str_num_proc);
  
  string output_dir = parser.getCmdOption("-o");
  if (output_dir.empty()){
    cout<<"missing output directory"<<endl;
    return 0;
  }
  
  //	int count = 10;
  //	int num_proc = 8;
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(num_proc);
  
  double t1 = omp_get_wtime();
  string line;
  
  cout<<"loading "<<indptr_file<<endl;
  ifstream findptr(indptr_file,ifstream::in);
  getline(findptr,line);
  vector<int> d_n;
  parse_string_2_int(line, d_n);
  int dim = d_n[0];
  int n = d_n[1];
  cout<<"n="<<n<<",dim="<<dim<<endl;
  
  vector<int> indptr;
  getline(findptr,line);
  parse_string_2_int(line, indptr);
  findptr.close();
  
  cout<<"loading "<<data_file<<endl;
  ifstream fdata(data_file,ifstream::in);
  getline(fdata,line);
  vector<float> data;
  parse_string_2_float(line, data);
  fdata.close();
  
  cout<<"loading "<<indices_file<<endl;
  ifstream findices(indices_file,ifstream::in);
  getline(findices,line);
  vector<int> indices;
  parse_string_2_int(line, indices);
  findices.close();
  t1 = omp_get_wtime() - t1;
  
  double t = omp_get_wtime();
  vector<float> dist(n, FLT_MAX); // will this cause a lot memory footprint? 
  vector<int> assign(n, -1);
  // 0 for non-landmark, 1 for landmark
  vector<int> flag(n, 0);
  select_landmarks_cpp(data, indices, indptr, n, dim, count, dist, assign, flag);
  t = omp_get_wtime() - t;
  
  // output files
  double t2 = omp_get_wtime();
  string landmark_file = output_dir+"/landmark.txt";
  ofstream out_landmark(landmark_file, ofstream::out);
  for(int i = 0;i<n;i++){
    if (flag[i] == 1){
      out_landmark<<i;
      for(int pos=indptr[i];pos<indptr[i+1];pos++){
        out_landmark<<" "<<indices[pos]<<":"<<data[pos];
      }
      out_landmark<<"\n";
    }
  }
  string assign_dist = output_dir+"/assign_dist.txt";
  ofstream out_dist(assign_dist, ofstream::out);
  for(int i=0;i<n;i++){
    out_dist<<i<<" "<<assign[i]<<" "<<dist[i]<<"\n";
  }
  t2 = omp_get_wtime() - t2;
  
  cout<<"loading data takes "<<t1<<" seconds"<<endl;
  cout<<"processing takes "<<t<<" seconds"<<endl;
  cout<<"output data takes "<<t2<<" seconds"<<endl;
  printf ("It took me %f seconds in total.\n",t1+t+t2);
  
  return 0;
}

// // [[Rcpp::export]]
// NumericMatrix jaccard_coeff(SEXP R_idx, SEXP R_weight) {
//   NumericMatrix idx(R_idx);
//   bool weight = as<bool>(R_weight);
//   
//   return jaccard_coeff_cpp(idx, weight);
// }

// [[Rcpp::export]]
Rcpp::List select_landmarks(SEXP R_data,
                            SEXP R_indices,
                            SEXP R_indptr,
                            SEXP R_n, // number of sample 
                            SEXP R_dim, // number of dimensioin 
                            SEXP R_count){ // number of landmarks
  
  std::vector<float> data = as<std::vector<float> >(R_data);
  std::vector<int> indices = as<std::vector<int> >(R_indices);
  std::vector<int> indptr = as<std::vector<int> >(R_indptr);
  int n = as<int>(R_n);
  int dim = as<int>(R_dim);
  int count = as<int>(R_count);
  vector<float> dist(n, FLT_MAX); // will this cause a lot memory footprint? 
  vector<int> assign(n, -1);
  // 0 for non-landmark, 1 for landmark
  vector<int> flag(n, 0);
  Rcpp::Rcout << "Before select_landmarks_cpp" << std::endl;
  
  // void select_landmarks_cpp(const vector<float> &data, const vector<int> &indices,
  //                           const vector<int> &indptr, int n, int dim, int count,
  //                           vector<float> &dist, vector<int> &assign, vector<int> &flag)
  
  select_landmarks_cpp(data, indices, indptr, n, dim, count, dist, assign, flag);
  Rcpp::Rcout << "after select_landmarks_cpp" << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("assign") = assign, //this really should be W. W can be used to for reverse embedding for missing data 
                            Rcpp::Named("dist") = dist,
                            Rcpp::Named("flag") = flag);
}

// test the code to identify the landmarks 
/*** R
library(monocle)
  lung <- load_lung()
  lung@assayData$exprs <- as(lung@assayData$exprs, 'sparseMatrix')  
  sp_data <- lung@assayData$exprs
  library(devtools)
  load_all()
  select_landmarks(sp_data@x, sp_data@i, sp_data@p, sp_data@Dim[2], sp_data@Dim[1], 10)
  */

/*** R
# library(monocle)
#   lung <- load_lung()

  
  library(devtools)
  load_all()
  
  lung@assayData$exprs <- as(lung@assayData$exprs, 'sparseMatrix')  
  sp_data <- lung@assayData$exprs
  system.time(select_landmarks(sp_data@x, sp_data@i, sp_data@p, sp_data@Dim[2], sp_data@Dim[1], 10))
  */

// test on Jun's dataset 
/*** R 
library(R.matlab)

load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle3/June/cds_mouse_eb_epi.RData') # 180201_mouse_embryo.RData 
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle3/June/180201_mouse_embryo.RData') # 180201_mouse_embryo.RData 
sp_data <- cds@assayData$exprs
start_time <- Sys.time()
Jun_landmark <- select_landmarks(sp_data@x, sp_data@i, sp_data@p, sp_data@Dim[2], sp_data@Dim[1], 1500)
end_time <- Sys.time()
end_time - start_time

Jun_X_1500 <- sp_data[, unique(Jun_landmark$assign) + 1] # converting into R index 
pData <- pData(cds)[unique(Jun_landmark$assign) + 1, ]
valid_Jun_X_1500 <- Jun_X_1500[, is.finite(pData[, 'Cluster'])]
valid_pData <- pData[is.finite(pData[, 'Cluster']), ]

writeMat('/Users/xqiu/Dropbox (Personal)/Projects/Monocle3/June/Jun_X_1500_mat.mat', Jun_X_1500 = valid_Jun_X_1500)

qplot(pData[, 'tsne_1'], pData[, 'tsne_2'], color = pData[, 'Cluster']) # reduce to top 30 pca, and then your method? 
qplot(pData[, 'tsne_1'], pData[, 'tsne_2'], color = pData[, 'day']) # reduce to top 30 pca, and then your method? 

write.table(pData(cds)[unique(Jun_landmark$assign) + 1, ], '/Users/xqiu/Dropbox (Personal)/Projects/Monocle3/June/Jun_X_1500_pd.csv', quote = F, sep = '\t')

irlba_res <- prcomp_irlba(t(log(valid_Jun_X_1500 + 1)), n = 20,
                          center = TRUE, scale. = TRUE)

irlba_res <- irlba((log(valid_Jun_X_1500 + 1)), nv = 20)
irlba_res <- pca_projection_R(log(valid_Jun_X_1500 + 1), L = 100)
writeMat('/Users/xqiu/Dropbox (Personal)/Projects/Monocle3/June/irlba_res.mat', Jun_X_1500 = t(irlba_res))

# run tSNE again: 
identical(row.names(pData(cds)), colnames(cds@assayData$exprs))
identical(row.names(fData(cds)), row.names(cds@assayData$exprs))
  
identical(row.names(pData(valid_cds)), colnames(valid_cds@assayData$exprs))
identical(row.names(fData(valid_cds)), row.names(valid_cds@assayData$exprs))

valid_cds <- cds[, !is.na(pData(valid_cds)[, 'tsne_2'])]
pd <- new("AnnotatedDataFrame", data = pData(valid_cds)[colnames(valid_cds@assayData$exprs), ])
fd <- new("AnnotatedDataFrame", data = fData(valid_cds))

valid_cds_corrected <- newCellDataSet(valid_cds@assayData$exprs, phenoData = pd, featureData = fd)

valid_cds <- estimateSizeFactors(valid_cds)
valid_cds <- estimateDispersions(valid_cds)

landmark_subset <- cds[, unique(colnames(valid_Jun_X_1500))]
landmark_subset <- estimateSizeFactors(landmark_subset)
landmark_subset <- estimateDispersions(landmark_subset)

Jun_skl <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/Monocle3/AAAI2017-code/matlab.mat')
qplot(Jun_skl$Y[, 1], Jun_skl$Y[, 2], color = pData(cds[, colnames(valid_Jun_X_1500)])$Cluster)

qplot(Jun_skl$Y[, 1], Jun_skl$Y[, 2], color = pData(cds[, unique(Jun_landmark$assign) + 1])$day) + xlab('SKL dim 1') + ylab('SKL dim 2') +  theme(legend.title = element_blank())
*/

/*** R
library(cellrangerRkit)
cellranger_pipestance_path <- "./"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)

fd <- fData(gbm)

# The number 2 is picked arbitrarily in the line below.
# Where "2" is placed you should place the column number that corresponds to your
# featureData's gene short names.

colnames(fd)[2] <- "gene_short_name"

gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
featureData = new("AnnotatedDataFrame", data = fd),
lowerDetectionLimit = 0.5,
expressionFamily = negbinomial.size())

 */