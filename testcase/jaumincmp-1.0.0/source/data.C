#include "data.h"
#include "math.h"

static bool CompareLevelData(LevelData& data1, LevelData& data2, double& tol,
                             Results& result) {
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
  for (size_t p = 0; p < data_2.size(); p++) {       // level
    for (size_t e = 0; e < data_2[p].size(); e++) {  // patch
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
            result.max_value1 = max(result.max_value1, v1);
            result.min_value1 = min(result.min_value1, v1);
            result.max_value2 = max(result.max_value2, v2);
            result.min_value2 = min(result.min_value2, v2);
            result.hits++;

            if (fabs(v1 - v2) > tol) {
              result.differ++;
              result.max_abs_error = max(fabs(v1 - v2), result.max_abs_error);
              result.min_abs_error = min(fabs(v1 - v2), result.min_abs_error);
              double relative_error = 0;
              if ((v1 + v2) != 0)
                relative_error = fabs(fabs(v1 - v2) / fabs(v1 + v2));
              result.max_rela_error =
                  max(result.max_rela_error, relative_error);
              result.min_rela_error =
                  min(result.max_rela_error, relative_error);
            } else {
              result.equal++;
            }
            pos++;
          }
        }
      }
      if (!found) {
        result.hits++;
        result.differ++;
      }
    }
  }

  return succ;
}

static bool CompareLevelDataTable(LevelDataTable& data1, LevelDataTable& data2,
                                  double& tol) {
  bool succ = true;

  Results result;
  result.hits = 0;
  result.equal = 0;
  result.differ = 0;
  result.max_abs_error = 0.0;
  result.min_abs_error = 10e300;
  result.max_rela_error = 0.0;
  result.min_rela_error = 10e300;
  result.max_value1 = -10e300;
  result.min_value1 = 10e300;
  result.max_value2 = -10e300;
  result.min_value2 = 10e300;

  // first check the variables of level data
  int nvar1 = data1.size();
  int nvar2 = data2.size();
  if (nvar1 != nvar2) {
    cerr << "VARIABLES_DIFFER.\n";
    return false;
  }
  // second compare each varibles data of levels
  map<string, LevelData>::iterator var_iter = data2.begin();
  while (var_iter != data2.end()) {
    if (!data1.count(var_iter->first)) {
      cerr << "VARIABLES_DIFFER.\n";
      return false;
    }
    succ =
        CompareLevelData(data1[var_iter->first], var_iter->second, tol, result);
    cout << var_iter->first << ": HITS=" << result.hits << ", ";
    if (result.differ == 0)
      cout << "TOTAL_EQUAL \n";
    else {
      cout << "EQUAL=" << result.equal << ", DIFFER=" << result.differ
           << ", MIN_ABS_ERROR=" << result.min_abs_error
           << ", MAX_ABS_ERROR=" << result.max_abs_error
           << ", MIN_RELATIVE_ERROR=" << result.min_rela_error
           << ", MAX_RELATIVE_ERROR=" << result.max_rela_error
           << ", MIN_VALUE1=" << result.min_value1
           << ", MAX_VALUE1=" << result.max_value1
           << ", MIN_VALUE2=" << result.min_value2
           << ", MAX_VALUE2=" << result.max_value2 << endl;
    }
    // if (!succ) {
    // 	cerr << "variable=" << var_iter->first << ", ";
    // 	return false;
    // }

    var_iter++;
  }

  return succ;
}

bool CompareHierarchyData(HierarchyData& data1, HierarchyData& data2,
                          double& tol) {
  bool succ = true;
  int nfed1 = data1.fed_data.size();
  int nfed2 = data2.fed_data.size();

  if (nfed1 != nfed2) {
    cerr << "FEDERAL_DIFFER_OR_LEVEL_DIFFER. \n";
    return false;
  }

  map<pair<int, int>, LevelDataTable>::iterator fed_iter =
      data2.fed_data.begin();
  while (fed_iter != data2.fed_data.end()) {
    cout << "-----------------------------------------\n";
    // if the data1 not contain the federation of data2
    if (!data1.fed_data.count(fed_iter->first)) {
      cerr << "FEDERAL_DIFFER_OR_LEVEL_DIFFER.\n";
      return false;
    } else {
      cout << "FEDERAL_EQUAL_ADD_LEVEL_EQUAL. \n";
    }

    LevelDataTable& level_table_data2 = fed_iter->second;
    LevelDataTable& level_table_data1 = data1.fed_data[fed_iter->first];

    succ = CompareLevelDataTable(level_table_data1, level_table_data2, tol);
    // if (!succ) {
    // 	cerr << "federal=" << fed_iter->first.first << endl;
    // 	return false;
    // }

    fed_iter++;
  }

  return succ;
}

void DeleteNullData(HierarchyData& data) {
  map<pair<int, int>, LevelDataTable>& federation_data = data.fed_data;
  map<pair<int, int>, LevelDataTable>::iterator fed_iter =
      federation_data.begin();
  while (fed_iter != federation_data.end()) {
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

void ExportData(HierarchyData& data, string filename) {
  ofstream fout(filename.c_str());
  map<pair<int, int>, LevelDataTable>& federation_data = data.fed_data;
  map<pair<int, int>, LevelDataTable>::iterator fed_iter =
      federation_data.begin();
  while (fed_iter != federation_data.end()) {
    map<string, LevelData>& levtabl = fed_iter->second;
    map<string, LevelData>::iterator var_iter = levtabl.begin();
    while (var_iter != levtabl.end()) {
      vector<vector<appu::CompareData<NDIM>::DataItem> >& level_data =
          var_iter->second.data;
      for (size_t p = 0; p < level_data.size(); p++) {
        vector<appu::CompareData<NDIM>::DataItem>& patch_data = level_data[p];
        for (size_t i = 0; i < patch_data.size(); i++) {
          int pos = 0;
          for (int d = 0; d < var_iter->second.meta_data.ndepth; d++)
            for (int g = 0; g < var_iter->second.meta_data.ndepth; g++) {
              fout << "federation:" << fed_iter->first.first
                   << " level:" << fed_iter->first.second
                   << " variable:" << var_iter->first << " patch:" << p
                   << " entity:" << i << "," << patch_data[i].entity_no
                   << " depth:" << d << " group:" << g;
              fout << "id=" << patch_data[i].id
                   << "cor1=" << patch_data[i].coord[0]
                   << "cor2=" << patch_data[i].coord[1];
              fout << " v=" << patch_data[i].data[pos] << endl;
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