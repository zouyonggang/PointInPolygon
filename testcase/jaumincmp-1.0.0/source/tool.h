#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "assert.h"
#include "hdf5.h"

using namespace std;


void GetChildFile(const string& path, vector<string>& files);
void GetChildGroup(const hid_t& group_id, vector<string>& child_group);

bool CompStepInSort(string a, string b);
void ReadDataset(const hid_t& group_id, const string& dataset_name, int& data);
void CheckDataset(const hid_t& group_id, const string& dataset_name) ;
double Dround(double number, unsigned int bits);
bool DirExist(string dirname);