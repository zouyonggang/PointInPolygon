#pragma once

#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <iomanip>
#include <fstream> 

#include "assert.h"
#include "hdf5.h"
#include "CompareData.h"
#include "PointOctree.h"
using namespace JAUMIN;
using namespace std;


//Stores data of a level
struct  LevelData
{
  appu::CompareData<NDIM>::Meta meta_data;
  vector<vector<appu::CompareData<NDIM>::DataItem> > data;
};

typedef map<string, LevelData> LevelDataTable;

// Stores data of all levels
struct HierarchyData
{
  map<pair<int,int>, LevelDataTable > fed_data;    
};

bool CompareHierarchyData(HierarchyData& data1, HierarchyData& data2, double& tol);
void ExportData(HierarchyData& data, string filename);
void DeleteNullData(HierarchyData& data);

#include "data.h"
#include "math.h"

static bool CompareLevelData(LevelData& data1, LevelData& data2, double& tol) 
{
    bool succ = true;
    appu::CompareData<NDIM>::Meta& meta_data_1 = data1.meta_data;
    appu::CompareData<NDIM>::Meta& meta_data_2 = data2.meta_data;
    vector<vector<appu::CompareData<NDIM>::DataItem> >& data_1 = data1.data;
    vector<vector<appu::CompareData<NDIM>::DataItem> >& data_2 = data2.data;

   // Check meta data
    assert(meta_data_1.data_id == meta_data_2.data_id);
    assert(meta_data_1.entity_type == meta_data_2.entity_type);
    assert(meta_data_1.ndepth == meta_data_2.ndepth);
    assert(meta_data_1.ngroup == meta_data_2.ngroup);

    // (1) Construct the octree with data_1
    // first I have to compute the bounding box
    double bbox[NDIM][2];
    for (int d = 0; d < NDIM; d++) {
      bbox[d][0] = 1.0e16;
      bbox[d][1] = -1.0e16;
    }
    for (size_t p = 0; p < data_1.size(); p++) {
      for (size_t e = 0; e < data_1[p].size(); e++) {
        for (int d = 0; d < NDIM; d++) {
          bbox[d][0] = std::min(bbox[d][0], data_1[p][e].coord[d]);
          bbox[d][1] = std::max(bbox[d][1], data_1[p][e].coord[d]);
        }
      }
    }
    // then new the tree
    tbox::PointOctree<NDIM, appu::CompareData<NDIM>::DataItem> tree(bbox);
    // insert the data to the tree
    for (size_t p = 0; p < data_1.size(); p++) {
      for (size_t e = 0; e < data_1[p].size(); e++) {
        tree.insert(data_1[p][e].coord, data_1[p][e]);
      }
    }

    // (2) Loop over data_2, compare against data_1
    int ndepth = meta_data_1.ndepth;
    int ngroup = meta_data_1.ngroup;
    for (size_t p = 0; p < data_2.size(); p++) {          //level
      for (size_t e = 0; e < data_2[p].size(); e++) {         //patch
        vector<appu::CompareData<NDIM>::DataItem> search_result;
        tree.search(data_2[p][e].coord, search_result);
        bool found = false;
        for (size_t t = 0; t < search_result.size(); t++) {
          if (search_result[t].id != data_2[p][e].id) continue;
          found = true;
          // check the data
          int pos = 0;
          for (int idepth = 0; idepth < ndepth; idepth++) {
            for (int igroup = 0; igroup < ngroup; igroup++) {
              double v1 = search_result[t].data[pos];
              double v2 = data_2[p][e].data[pos];
              if (fabs(v1 - v2) > tol) {
                // fprintf(stderr,
                //         "fail at comparing <patch=%d,entity=%d> of dir2,  with "
                //         "<patch=%d,entity=%d> of dir1\n",
                //         p, e, search_result[t].patch_no,
                //         search_result[t].entity_no);
                cerr << "fail at comparing <patch=" << p << ",entity=" << e
                         << "> of dir2,  with <patch=" << search_result[t].patch_no 
                         << ",entity=" << search_result[t].entity_no << "> of dir1\n";
                fprintf(stderr,
                        "v1=%0.20le, v2=%0.20le, err=%0.4le, tol=%0.4le, ", v1,
                        v2, fabs(v1 - v2), tol);
                return false;
              }
              pos++;
            }
          }
        }
        if (!found) {
          // fprintf(stderr, "can't find <patch=%d,entity=%d> of dir2 in dir 1.\n",
          //         p, e);
          cerr << "can't find <patch=" << p << ",entity=" << e <<"> of dir2 in dir 1.\n";
          return false;
        }
      }
    }

return succ;
}

static bool CompareLevelDataTable(LevelDataTable& data1, LevelDataTable& data2, double& tol) 
{
  bool succ = true;

  //first check the variables of level data
  int nvar1 = data1.size();
  int nvar2 = data2.size();
  if (nvar1 != nvar2) {
    cerr << "mismatch the number of variables between levels.\n";
    return false;
  }
  //second compare each varibles data of levels
  map<string, LevelData>::iterator var_iter = data2.begin();
  while (var_iter != data2.end()) {
    if ( !data1.count(var_iter->first) ) {
      cerr << "mismatch the number of variables between levels\n";
      return false;
    }
    succ = CompareLevelData(data1[var_iter->first], var_iter->second, tol);
    if (!succ) {
      cerr << "variable=" << var_iter->first << ", ";
      return false;
    }

    var_iter++;
  }
  
  return succ;
}

bool CompareHierarchyData(HierarchyData& data1, HierarchyData& data2, double& tol)
{
  bool succ = true;
  int nfed1 = data1.fed_data.size();
  int nfed2 = data2.fed_data.size();
  
  if (nfed1 != nfed2) {
    cerr << "mismatch the number of federation. \n";
    return false;
  }

  map<pair<int,int>, LevelDataTable >::iterator fed_iter = data2.fed_data.begin();
  while (fed_iter != data2.fed_data.end()) {
    //if the data1 not contain the federation of data2
    if (! data1.fed_data.count(fed_iter->first)) {
      cerr << "mismatch the federation of data1 and data2.\n";
      return false;
    } 

    LevelDataTable& level_table_data2 = fed_iter->second;
    LevelDataTable& level_table_data1 = data1.fed_data[fed_iter->first];

    succ = CompareLevelDataTable(level_table_data1, level_table_data2, tol);
    if (!succ) {
      cerr << "federal=" << fed_iter->first.first << endl;
      return false;
    }

    fed_iter ++;
  } 

  return succ;
}

void DeleteNullData(HierarchyData& data)
{
  map<pair<int,int>, LevelDataTable >& federation_data = data.fed_data;
  map<pair<int,int>, LevelDataTable >::iterator fed_iter = federation_data.begin();
  while(fed_iter != federation_data.end()) {
    if ((fed_iter->second).size() == 0)
      federation_data.erase(fed_iter++);
    else {
      map<string, LevelData>& levtabl = fed_iter->second;
      map<string, LevelData>::iterator var_iter = levtabl.begin();
      while (var_iter != levtabl.end()) {
        if ((var_iter->second).data.size() == 0)
          levtabl.erase(var_iter++);
        else 
          var_iter++;
      }
      fed_iter++;
    }
  }
}

void ExportData(HierarchyData& data, string filename)
{
  ofstream fout(filename.c_str());
  map<pair<int,int>, LevelDataTable >& federation_data = data.fed_data;
  map<pair<int,int>, LevelDataTable >::iterator fed_iter = federation_data.begin();
  while(fed_iter != federation_data.end()) {
    map<string, LevelData>& levtabl = fed_iter->second;
    map<string, LevelData>::iterator var_iter = levtabl.begin();
    while (var_iter != levtabl.end()) {
      vector<vector<appu::CompareData<NDIM>::DataItem> >& level_data
        = var_iter->second.data;
        for (size_t p = 0; p < level_data.size(); p++) {
          vector<appu::CompareData<NDIM>::DataItem>& patch_data = level_data[p];
          for (size_t i = 0; i < patch_data.size(); i++) {
            int pos = 0;
            for(int d = 0; d < var_iter->second.meta_data.ndepth; d++)
              for(int g = 0; g < var_iter->second.meta_data.ndepth; g++){
                fout << "federation:" << fed_iter->first.first << " level:" 
                  << fed_iter->first.second << " variable:" << var_iter->first 
                  << " patch:" << p << " entity:" << i <<","<< patch_data[i].entity_no 
                  << " depth:" << d << " group:" << g;
                fout <<setprecision(20)<<"id="<<patch_data[i].id<< "cor1=" << patch_data[i].coord[0]<<"cor2="<<patch_data[i].coord[1];
                fout << setprecision(20) << " v=" << patch_data[i].data[pos] << endl;
                pos++;
              }
          } 
        }
       var_iter++;
    }
    fed_iter++;
  }
  fout.close();
}
















// This code compares the outputs (either in javis format
// or restart format) of different runs of a JAUMIN application.
//
// Usage:
//     validation -dir1 <dirname> -dir2 <dirname> [other_options]
// other_options can be:
//     -format <javis|restart>
//     -tol <tolerance>
//     -verbose <0|1>
//     -only_last_step
//     -export_data
//
// Author: Zou Yonggang, Cheng Jie, 2017.


#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <stdexcept>

#include "JAUMIN_config.h"
#include "JAUMINManager.h"
#include "MPI.h"

#include "restart.h"
#include "javis.h"
#include "tool.h"
#include "hdf5.h"

using namespace std;
using namespace JAUMIN;

struct Params {
  string datadir1;
  string datadir2;
  double tol;
  int verbose;
  string format;  // "javis" | "restart"
  bool only_last_step;
  bool export_data;
  string app_name;
  string config_file;
} params; /* 存储数据比对的初始数据 */

void ParseParams(int argc, char* argv[]) {
  // Set the default value
  params.verbose = 1;
  params.format = "restart";
  params.only_last_step = false;
  params.tol = 1.0e-15;
  params.app_name = "unknown";
  params.export_data = false;

  int pos = 1;
  while (pos < argc) {
    string key(argv[pos]);
    if (key == "-dir1") {
      params.datadir1 = argv[pos + 1];
      pos += 2;
    } else if (key == "-dir2") {
      params.datadir2 = argv[pos + 1];
      pos += 2;
    }  else if (key == "-verbose") {
      params.verbose = atoi(argv[pos + 1]);
      pos += 2;
    } else if (key == "-format") {
      params.format = argv[pos + 1];
      pos += 2;
    } else if (key == "-only_last_step") {
      params.only_last_step = true;
      pos += 1;
    } else if (key == "-export_data") {
      params.export_data = true;
      pos += 1;
    }
  else if(key == "-tol"){
    stringstream ss(argv[pos + 1]);
    ss >> params.tol;
    pos += 2;
    }
  else if(key == "-app_name"){
    params.app_name = argv[pos + 1];
    pos += 2;
    } 
  else {
     cerr << "please input such as :" << endl;
      cerr << "   ./main2d -dir1 /path/to/src/restart/dir -dir2 /path/to/target/restart/dir" << endl;
      throw invalid_argument("unknown key word.\n");
    }
  }
  cout << "dir1: " << params.datadir1 << endl;
  cout << "dir2: " << params.datadir2 << endl;
  cout << "tol : " << params.tol << endl;
}

bool CompareData(string& data_dir_1, string& data_dir_2) {
  vector<string> stepdirs1, stepdirs2;

  if (params.format == "restart") {
    CollectRestartStepDir(data_dir_1, stepdirs1);
    CollectRestartStepDir(data_dir_2, stepdirs2);
  } else if (params.format == "javis") {
    CollectJavisStepDir(data_dir_1, stepdirs1);
    CollectJavisStepDir(data_dir_2, stepdirs2);
  }
  sort(stepdirs1.begin(),stepdirs1.end(),CompStepInSort);
  sort(stepdirs2.begin(),stepdirs2.end(),CompStepInSort);

  // First, check whther the step dirs are the same
  if (stepdirs1.size() != stepdirs2.size()) {
    if (params.verbose > 0) cerr << "mismatch number of time steps\n";
    return false;
  }

  // If there are no data, consider it succ.
  if (stepdirs1.size() == 0u) return true;

  if (params.only_last_step) {
    int n = (int)stepdirs1.size();
    if (stepdirs1[n - 1] != stepdirs2[n - 1]) {
      if (params.verbose > 0) cerr << "mismatch of time step\n";
      return false;
    }
  stepdirs1[n-1] = data_dir_1 + "/" + stepdirs1[n-1];
  stepdirs2[n-1] = data_dir_2 + "/" + stepdirs2[n-1];
  } else {
    int n = (int)stepdirs1.size();
    for (int i = 0; i < n; i++) {
      if (stepdirs1[i] != stepdirs2[i]) {
        if (params.verbose > 0) {
          cerr << "mismatch of time steps\n";
          cerr << stepdirs1[i] << " v.s. " << stepdirs2[i] << endl;
        }
        return false;
      }
    stepdirs1[i] = data_dir_1 + "/" + stepdirs1[i];
    stepdirs2[i] = data_dir_2 + "/" + stepdirs2[i];
    }
  }

  // Second, check the data of each step
  
  // Set the compare params
  //CompareParam compare_parameters;
  //ReadCompareParameters(config_file, app_name, compare_parameters)
  
  bool succ = true;
  
  if (params.only_last_step) {
    int n = (int)stepdirs1.size();
    if (params.format == "restart") {
       succ = CompareRestartStepData(stepdirs1[n - 1], stepdirs2[n - 1], params.tol);
    } else if (params.format == "javis") {
      succ = CompareJavisStepData(stepdirs1[n - 1], stepdirs2[n - 1]);
    }
  } else {
    int n = (int)stepdirs1.size();
  if (params.format == "restart") {
    for (int i = 0; i < n; i++) {
      succ = CompareRestartStepData(stepdirs1[i], stepdirs2[i], params.tol);
      if (!succ) break;
    }
  } else if (params.format == "javis") { //the first file is the summary
    for (int i = 1; i < n; i++) {
      succ = CompareJavisStepData(stepdirs1[i], stepdirs2[i]);
      if (!succ) break;
    }
  }
  }

  return succ;
}

int main(int argc, char* argv[]) {
  tbox::MPI::init(&argc, &argv);
  tbox::JAUMINManager::startup();
  ParseParams(argc, argv);
  bool succ = CompareData(params.datadir1, params.datadir2);
  if (succ) 
  cout<<"succ=true"<<endl;
  else
    cout<<"succ=false"<<endl;
  // 释放JAUMIN和MPI内部资源.
  tbox::JAUMINManager::shutdown();
  tbox::MPI::finalize();
  return (succ ? 0 : 1);
}
