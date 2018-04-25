#pragma once

#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <iomanip>
#include <fstream> 
#include <algorithm>

#include "assert.h"
#include "hdf5.h"
#include "CompareData.h"
#include "PointOctree.h"
using namespace JAUMIN;
using namespace std;

struct  Results
{
	int hits;
	int equal;
	int differ;
	double max_abs_error;
	double min_abs_error;
	double max_rela_error;
	double min_rela_error;
	double max_value1;
	double min_value1;
	double max_value2;
	double min_value2;
};

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