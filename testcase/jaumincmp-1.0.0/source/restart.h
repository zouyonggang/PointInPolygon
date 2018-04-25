#pragma once

#include <string>
#include <vector>

using namespace std;


void CollectRestartStepDir(const string& data_dir_1, vector<string>& stepdirs1);
bool CompareRestartStepData(string& step_dir_1, string& step_dir_2, double& tol );


