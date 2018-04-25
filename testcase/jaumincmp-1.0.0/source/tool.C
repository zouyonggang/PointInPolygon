#include <dirent.h>
#include "tool.h"
#include <cstring>
#include <sstream>
#include <iomanip>
#include "sys/stat.h"

void GetChildFile(const string& path, vector<string>& files)
{
    DIR              *pDir ; 
    struct dirent    *ent  ;
  
    pDir=opendir(path.c_str());
  
    while((ent=readdir(pDir))!=NULL)
    {  
		if(strcmp(ent->d_name,".")==0 || strcmp(ent->d_name,"..")==0)
			continue;
		string temp;
		temp = ent->d_name;
		//temp = path + '/' + temp;
		files.push_back(temp);
    }	
}

void GetChildGroup(const hid_t& group_id, vector<string>& child_group)
{
	int i,MAX_NAME=64;
	ssize_t len;
	hsize_t nobj;
	herr_t err;
	int otype;
	char child_group_name[MAX_NAME];
	
	err = H5Gget_num_objs(group_id,&nobj);
	if((int)err < 0) cerr << "open hdf5 dataset fail.\n";
	
	for(i = 0; i < (int)nobj; i++){
		len = H5Gget_objname_by_idx(group_id, (hsize_t)i, 
		child_group_name, (size_t)MAX_NAME );
		if((int)len < 0) std::cout << "open hdf5 group fail.\n";
		
		otype = H5Gget_objtype_by_idx(group_id,(size_t)i);
		
		if(otype == H5G_GROUP){
			child_group.push_back(child_group_name);
		}
	}
}

//check the dataset exit
void CheckDataset(const hid_t& group_id, const string& dataset_name) 
{
  htri_t exit_status = H5Lexists(group_id, dataset_name.c_str(), H5P_DEFAULT);
  assert(exit_status >= 0);
}

//ReadDataset
void ReadDataset(const hid_t& group_id, const string& dataset_name, int& data)
{
	hid_t dataset;
	herr_t status;
	CheckDataset(group_id, dataset_name);
	dataset = H5Dopen2(group_id, dataset_name.c_str(), H5P_DEFAULT);
	assert(dataset > 0);
 	status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&data);
	assert(status >= 0);
 	H5Dclose(dataset);
}

//sort the vector which store the stepdirs
bool CompStepInSort(string a, string b) 
{
	int a1 = atoi(a.substr(8,5).c_str());
	int b1 = atoi(b.substr(8,5).c_str());
	return a1 < b1;
}

double Dround(double number, unsigned int bits) {
    stringstream ss;
    ss << fixed << setprecision(bits) << number;
    ss >> number;
    return number;
}

bool DirExist(string dir)
{
	struct stat fileStat;
	if ((stat(dir.c_str(), &fileStat) == 0) && S_ISDIR(fileStat.st_mode))
		return true;
	else {
		cerr << dir << "is not exist!" <<endl;
		return false;
	}
}