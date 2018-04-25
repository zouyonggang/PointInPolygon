#pragma once

#include <string>
#include <vector>

using namespace std;


void CollectJavisStepDir(const string& data_dir_1, vector<string>& stepdirs1);
bool CompareJavisStepData(string& step_dir_1, string& step_dir_2);