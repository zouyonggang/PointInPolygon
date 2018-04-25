#include "restart.h"
#include "data.h"
#include "tool.h"

// ---------------- private utilities----------------

//Get variables from level
static void GetVariables(const hid_t& patch_group_id, vector<string>& variable_list)
{
  hid_t dataset, datatype, dataspace;
  size_t datasize;
  herr_t status;
  hsize_t number_select;

  htri_t exit_status = H5Lexists(patch_group_id, "patch_data_namelist", H5P_DEFAULT);
  if (!(exit_status >= 0) )
    return ;
  dataset = H5Dopen2(patch_group_id, "patch_data_namelist", H5P_DEFAULT);
  datatype = H5Dget_type(dataset);
  datasize = H5Tget_size(datatype);       //获取数据集中元素的大小
  dataspace = H5Dget_space(dataset);
  number_select = H5Sget_select_npoints(dataspace);//获取数据集中元素（字符）的个数
  
  //从hdf数据集中读取数据
  char* src_data;
  src_data = new char[datasize*number_select];
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL,H5P_DEFAULT,src_data);
  if ((int)status < 0) cerr<<"open hdf5 dataset fail.\n";
  
  //整理读取到的数据
  variable_list.resize((int)number_select);
  for (int i = 0; i < (int)number_select; i++) {
    string *locPtr = &variable_list[i];
    *locPtr = &src_data[i*datasize];
  }
  for(vector<string>::iterator it = variable_list.begin(); it != variable_list.end();) {
    if( (*it).find("label") == string::npos &&
          (*it).find("id") == string::npos )
      it = it + 1;
    else
      variable_list.erase(it);
  }

  H5Dclose(dataset);
  H5Tclose(datatype);
  H5Sclose(dataspace);

}

//Get meta_data from level's group
static void GetMetaData(const hid_t& level_group_id,  
                              appu::CompareData<NDIM>::Meta& meta_data,
                              const string& variable,
                              const string& patch_name)
{
  meta_data.data_id = 1;
  meta_data.ghost_width = 1;
  //get entity type
  ReadDataset(level_group_id, patch_name + "/" + variable + "/d_type", meta_data.entity_type);
  //get the number of patch
  ReadDataset(level_group_id, "d_number_global_patches", meta_data.npatch);
  //get ndepth
  ReadDataset(level_group_id, patch_name + "/" + variable + "/d_data/d_depth", meta_data.ndepth);
  //get ngroup
  ReadDataset(level_group_id, patch_name + "/" + variable + "/d_data/d_group", meta_data.ngroup);
}

//get cell geometry from adjacency_nodes_extent and adjacency_nodes_index
static void GetPhysicsFieldGeometry(const hid_t& patch_group_id, 
                                  string fieldname,
                                  vector<vector<double> >& coords, 
                                  vector<string>& item_ids)
{
  string dataset_name = fieldname + "/d_adjacency_nodes_extent";
  CheckDataset(patch_group_id, dataset_name);

//first:get adjacency_nodes_extent and adjacency_nodes_index
  hid_t dataset, dataspace;
  herr_t status;
  hsize_t extent_number,index_number;
  //get nodes extenn

  dataset = H5Dopen2(patch_group_id, dataset_name.c_str(), H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  extent_number = H5Sget_select_npoints(dataspace);
  int nodes_extent[extent_number];
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nodes_extent);
  assert(status >= 0);
  //get nodes index
  dataset_name = fieldname + "/d_adjacency_nodes_index";
  dataset = H5Dopen2(patch_group_id, dataset_name.c_str(), H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  index_number = H5Sget_select_npoints(dataspace);
  int nodes_index[index_number];
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nodes_index);
  assert(status >= 0);

//second:get cell Geometry
vector<vector<double> > cell_coords;
vector<string> cell_ids;
for (int e = 0; e < (int)extent_number-1; e++) {
  vector<double> cell_coord(3, 0);
  char cell_id[8] = {'z','z','z','z','z','z','z','z'};
  for (int i = nodes_extent[e]; i < nodes_extent[e+1]; i++) {
    for (int dim = 0; dim < NDIM; dim++) 
      cell_coord[dim] += coords[nodes_index[i]][dim];
    for (int c = 0; c < 8; c++)
      if (cell_id[c] > item_ids[nodes_index[i]][c])
        cell_id[c] = item_ids[nodes_index[i]][c];
  }
  for (int dim = 0; dim < NDIM; dim++) {
    cell_coord[dim] = cell_coord[dim]/(nodes_extent[e+1] - nodes_extent[e]);
    cell_coord[dim] = Dround(cell_coord[dim], 15);
  }
  cell_coords.push_back(cell_coord);
  string temp(cell_id, 8);
  cell_ids.push_back(temp);
}
coords.swap(cell_coords);
item_ids.swap(cell_ids);
vector<vector<double> >().swap(cell_coords);
vector<string>().swap(cell_ids);
}

//Get patch Geometry 
static void GetPatchGeometry(const hid_t& patch_group_id, int& entity_type, 
                                    vector<vector<double> >& coords,
                                    vector<string>& item_ids)
{
  hid_t dataset, datatype;
  size_t datasize;
  herr_t status;

 //first: get nodes coords and item_ids
  //get coordinates
  CheckDataset(patch_group_id, "d_patch_geometry/d_node_coordinates/d_data/d_array");
  int data_num;
  ReadDataset(patch_group_id, "d_patch_geometry/d_node_coordinates/d_data/d_entity", data_num);
  dataset = H5Dopen2(patch_group_id, 
      "d_patch_geometry/d_node_coordinates/d_data/d_array", H5P_DEFAULT);
  double coordinates[data_num * NDIM];
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &coordinates);
  assert(status >= 0);
  H5Dclose(dataset);

  //get ids
  CheckDataset(patch_group_id, "d_patch_geometry/d_node_ids/d_data/d_array");
  dataset = H5Dopen2(patch_group_id, "d_patch_geometry/d_node_ids/d_data/d_array", H5P_DEFAULT);
  datatype = H5Dget_type(dataset);
  datasize = H5Tget_size(datatype);
  char* ids;
  ids = new char[datasize];
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
  assert(status >= 0);
  H5Dclose(dataset);
  H5Tclose(datatype);

//second: deal with geogetry
  for (int i = 0; i < data_num*NDIM; i++) 
    coordinates[i] = Dround(coordinates[i], 20);

  item_ids.resize(data_num);
  for (int i = 0; i < data_num; i++) {
    char str[8];
    strncpy(str, ids+i*8, 8);
    string tmp(str,8);
    item_ids[i] = tmp;

   vector<double> temp;
#if (NDIM == 2)
    temp.push_back(coordinates[i*2]);
    temp.push_back(coordinates[i*2+1]);
    coords.push_back(temp);
#else
    temp.push_back(coordinates[i*3]);
    temp.push_back(coordinates[i*3+1]);
    temp.push_back(coordinates[i*3+2]);
    coords.push_back(temp);
#endif
    vector<double>().swap(temp);
  }

  if (entity_type == 1) //it's node data
    return;
  else if (entity_type == 0) { //it's cell data
    GetPhysicsFieldGeometry(patch_group_id, "d_patch_topology/d_cell_adjacency", coords, item_ids);
  }
  else if (entity_type == 2) { //it's edge data
    GetPhysicsFieldGeometry(patch_group_id, "d_patch_topology/d_edge_adjacency", coords, item_ids);
  }
  else if (entity_type == 3) { //it's face data
    GetPhysicsFieldGeometry(patch_group_id, "d_patch_topology/d_face_adjacency", coords, item_ids);
  }
  else
    cerr << "unknow data type, deal with like nodes" << endl;

}

//Read patch data from patch group
static void ReadPatchData(const hid_t& patch_group_id, 
                                vector<appu::CompareData<NDIM>::DataItem>& patch_data,
                                string& variable,
                                int& entity_type, int& depth, int& group) 
{
  hid_t dataset;
  herr_t status;
  string dataset_name;
  int d_patch_index, d_entity, d_offset;

//first:get Geometry
  vector<vector<double> > coords;
  vector<string> item_ids;
  GetPatchGeometry(patch_group_id, entity_type, coords, item_ids);
  assert(coords.size() == item_ids.size());

//second:get data
  //get patch's index
  ReadDataset(patch_group_id, "d_index", d_patch_index);

  //get entity
  ReadDataset(patch_group_id, variable + "/d_data/d_entity", d_entity);
  assert(d_entity <= (int)item_ids.size());

  //get offset
  ReadDataset(patch_group_id, variable + "/d_data/d_offset", d_offset);
  //assert(d_offset == d_entity*depth*group);

  //get data array
  dataset_name = variable + "/d_data/d_array";
  dataset = H5Dopen2(patch_group_id, dataset_name.c_str(), H5P_DEFAULT);
  assert(dataset > 0);
  hid_t dataspace = H5Dget_space(dataset);
  hsize_t number_select = H5Sget_select_npoints(dataspace);//获取数据集中元素（字符）的个数
  assert((int)number_select == d_entity*depth*group);
  double d_data[d_entity*depth*group];
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &d_data);
  assert(status >= 0);

  //deal with data
  for (int i = 0; i < d_entity; i++) {
    appu::CompareData<NDIM>::DataItem item_data;
    for (int dim = 0; dim < NDIM; dim++)
      item_data.coord[dim] = coords[i][dim];
    item_data.id = item_ids[i];
    item_data.patch_no = d_patch_index;
    item_data.entity_no = i;
    for (int d = 0; d < depth; d++) {
      int pos = i*group;
      for (int g = 0; g < group; g++) {
        item_data.data.push_back(d_data[pos + d_entity*group*d]);
        pos ++;
      }
    }
    patch_data.push_back(item_data);
  }
  vector<vector<double> >().swap(coords);
  vector<string>().swap(item_ids);
}

//Read Level data from level's group in proc file
static void ReadLevelData(const hid_t& level_group_id, LevelDataTable& levtab)
{
  //get patch list
  vector<string> patch_list;
  GetChildGroup(level_group_id, patch_list);
  for (vector<string>::iterator i = patch_list.begin(); i != patch_list.end();)
    if ((*i).find("patch") == string::npos)
      patch_list.erase(i);
    else
      i++;
  if (patch_list.size() == 0) return;

  //get varibles from level
  hid_t patch_group_id = 
    H5Gopen2(level_group_id, patch_list[0].c_str(), H5P_DEFAULT);
  vector<string> variable_list;
  GetVariables(patch_group_id, variable_list);
  H5Gclose(patch_group_id);
  if ((int)variable_list.size() <= 0) return;

  for (size_t v = 0; v < variable_list.size(); v++) {
    LevelData level_data;
    //get meta_data of level's
    GetMetaData(level_group_id, level_data.meta_data, variable_list[v], patch_list[0]);

    // if (level_data.meta_data.entity_type == 2 || 
    //       level_data.meta_data.entity_type == 3)
    //   continue;

    for (size_t p = 0; p < patch_list.size(); p++) {
      vector<appu::CompareData<NDIM>::DataItem> patch_data;
      patch_group_id = 
        H5Gopen2(level_group_id, patch_list[p].c_str(), H5P_DEFAULT);
      assert(patch_group_id > 0);
      ReadPatchData(patch_group_id, patch_data, variable_list[v],
                        level_data.meta_data.entity_type,
                        level_data.meta_data.ndepth,
                        level_data.meta_data.ngroup);
      H5Gclose(patch_group_id);
      if (patch_data.size() > 0)
        level_data.data.push_back(patch_data);
      vector<appu::CompareData<NDIM>::DataItem>().swap(patch_data);
    }
    levtab[variable_list[v]] = level_data;
  }
}

//Read hierarcy data from Hierarchy group in proc file
static void ReadHierarchyData(const hid_t& hierarchy_group_id, HierarchyData& data)
{
  //get child group of hierarchy group
  vector<string> hierarchy_child_list;
  GetChildGroup(hierarchy_group_id, hierarchy_child_list);
  assert(hierarchy_child_list.size() > 0);

  // if the child group is level
  if (hierarchy_child_list[0].find("level") != string::npos) {
    hid_t level_group_id;
    for (size_t l = 0; l < hierarchy_child_list.size(); l++){
      LevelDataTable levtable;
      level_group_id = 
        H5Gopen2(hierarchy_group_id, hierarchy_child_list[l].c_str(), H5P_DEFAULT);
      assert(level_group_id > 0);
      ReadLevelData(level_group_id, levtable);
      if (levtable.size() > 0)
        data.fed_data[make_pair(0,l)] = levtable;
      H5Gclose(level_group_id);
    }
  }
  //the child group is federal
  else if (hierarchy_child_list[0].find("federal") != string::npos) {
    hid_t federation_group_id;
    for (size_t f = 0; f <hierarchy_child_list.size(); f++) {
      federation_group_id = 
        H5Gopen2(hierarchy_group_id, hierarchy_child_list[f].c_str(), H5P_DEFAULT);
      vector<string> federation_child_list;
      GetChildGroup(federation_group_id, federation_child_list);
      assert(federation_child_list.size() > 0);

      hid_t level_group_id;
      for (size_t l = 0; l < federation_child_list.size(); l++) {
        level_group_id = 
          H5Gopen2(federation_group_id, federation_child_list[l].c_str(), H5P_DEFAULT);
        assert(level_group_id > 0);
        LevelDataTable levtable;
        ReadLevelData(level_group_id, levtable);
        if (levtable.size() > 0) 
          data.fed_data[make_pair(f,l)] = levtable;
        H5Gclose(level_group_id);
      }

      H5Gclose(federation_group_id);
    }
  } else {
      cerr << "unknow datatype: level? federal?" << endl;
      return;
    }

}

//Read a process data from proc file
static void ReadProcData(const string& file_name, HierarchyData& data)
{
  //读取每一个进程文件数据
  hid_t file,patch_hier_group_id;
  file = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file > 0);

  //找到进程文件根群组下存储数据的PatchHierarchy群组
  vector<string> root_child_list;
  GetChildGroup(file, root_child_list);
  for (int i = 0; i < (int)root_child_list.size(); i++)
    if (root_child_list[i].find("PatchHierarchy") != string::npos) {
      //打开PatchHierarchy读取数据
      patch_hier_group_id = H5Gopen2(file, root_child_list[i].c_str(), H5P_DEFAULT);
      assert(patch_hier_group_id > 0);
      ReadHierarchyData(patch_hier_group_id, data);
      H5Gclose(patch_hier_group_id);
    }
    
  H5Fclose(file);
}

//Merge the patches
static void MergeProcData(HierarchyData& root_data, HierarchyData& data)
{
  map<pair<int,int>, LevelDataTable >& fed_data1 = root_data.fed_data;
  map<pair<int, int>, LevelDataTable >& fed_data2 = data.fed_data;

  map<pair<int,int>, LevelDataTable >::iterator federal_iter = fed_data2.begin();
  while (federal_iter != fed_data2.end() ) {
    if (fed_data1.count(federal_iter->first) == 0)
      fed_data1[federal_iter->first] = federal_iter->second;
    else {
      map<string, LevelData>& levtab1 = fed_data1[federal_iter->first];
      map<string, LevelData>& levtab2 = federal_iter->second;

      map<string, LevelData>::iterator var_iter = levtab2.begin();
      while (var_iter != levtab2.end() ) {
        if (levtab1.count(var_iter->first) ==0 )
          levtab1[var_iter->first] = var_iter->second;
        else {
          LevelData& level_data1 = levtab1[var_iter->first];
          LevelData& level_data2 = var_iter->second;

          assert(level_data1.meta_data.entity_type == level_data1.meta_data.entity_type);
          assert(level_data1.meta_data.npatch == level_data1.meta_data.npatch);
          assert(level_data1.meta_data.ndepth == level_data1.meta_data.ndepth);
          assert(level_data1.meta_data.ngroup == level_data1.meta_data.ngroup);

          for (size_t p = 0; p < level_data2.data.size(); p++) 
            level_data1.data.push_back(level_data2.data[p]);
        }
        var_iter++;
      }
    }
  federal_iter++;
  }
}

// Read everything from the restart file
// Store the data in 'data' (see data.h)
static void ReadRestartData(const string& step_dir, HierarchyData& data)
{
  vector<string> node_list;
  GetChildFile(step_dir, node_list);
  node_list[0] = step_dir + "/" + node_list[0];

  //获取该节点下的所有进程文件
  vector<string> proc_list;
  GetChildFile(node_list[0], proc_list);

  //遍历所有的进程文件获取每个进程的数据
  for (int i = 0; i < (int)proc_list.size(); i++) {
      if (proc_list[i].find("proc") == string::npos)
      continue;
    proc_list[i] = node_list[0] + "/" + proc_list[i];
    if (i == 0) {//如果是第一个proc,则将数据取出来直接放进data
      ReadProcData(proc_list[i], data);
    } else {//如果不是第一个proc，则将数据取出来插入data,合并数据
      HierarchyData temp;
      ReadProcData(proc_list[i], temp);
      MergeProcData(data, temp);
    }
  }

}

// ---------------- public interface ----------------

void CollectRestartStepDir(const string& data_dir_1, vector<string>& stepdirs1)
{
  // collect the directories like restore.00010
	GetChildFile(data_dir_1, stepdirs1);
}

bool CompareRestartStepData(string& step_dir_1, string& step_dir_2, double& tol)
{
  HierarchyData d1, d2;
  //cout << "Read data1 ...\n";
  ReadRestartData(step_dir_1, d1);
  DeleteNullData(d1);
  //ExportData(d1, "data1.txt");
  //cout << "Read data2 ...\n";
  ReadRestartData(step_dir_2, d2);
  DeleteNullData(d2);
  //ExportData(d2, "data2.txt");
  //cout << "Compare data1 and data2 ...\n";
  bool r = CompareHierarchyData(d1, d2, tol);
  return r;
}
