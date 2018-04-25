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

  if (argc < 4) {
      cerr << "please input such as :" << endl;
      cerr << "   ./main2d /path/to/src/restart_dir/restore001 /path/to/target/restart_dir/restore001 1.0e-10" << endl;
      throw invalid_argument("unknown key word.\n");
  }
  params.datadir1 = argv[1];
  params.datadir2 = argv[2];
  stringstream ss(argv[3]);
  ss >> params.tol;

  cout << "dir1: " << params.datadir1 << endl;
  cout << "dir2: " << params.datadir2 << endl;
  cout << "tol : " << params.tol << endl;
}

bool CompareData(string& data_dir_1, string& data_dir_2) {
  if (! DirExist(data_dir_1) )
    return false;
  if (! DirExist(data_dir_2) )
    return false;
  return CompareRestartStepData(data_dir_1, data_dir_2, params.tol);
}

int main(int argc, char* argv[]) {
  tbox::MPI::init(&argc, &argv);
  tbox::JAUMINManager::startup();
  ParseParams(argc, argv);
  bool succ = CompareData(params.datadir1, params.datadir2);

  // 释放JAUMIN和MPI内部资源.
  tbox::JAUMINManager::shutdown();
  tbox::MPI::finalize();
  return (succ ? 0 : 1);
}
