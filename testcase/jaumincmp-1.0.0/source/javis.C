#include "tool.h"
#include "javis.h"
#include "data.h"

// ---------------- private utilities----------------

static void read_hierarchy_data(const string& dir, HierarchyData& data)
{
  // Read everything from the restart file
  // Store the data in 'data' (see data.h)
  //data.level_data.resize(0);
}

// ---------------- public interface ----------------

void CollectJavisStepDir(const string& data_dir_1, vector<string>& stepdirs1)
{
  // collect the directories like "javis_dump.00010"
  GetChildFile(data_dir_1,stepdirs1);
}

bool CompareJavisStepData(string& step_dir_1, string& step_dir_2)
{
  HierarchyData d1, d2;
  cout << "Read data1 ...\n";
  read_hierarchy_data(step_dir_1, d1);
  cout << "Read data2 ...\n";
  read_hierarchy_data(step_dir_2, d2);
  cout << "Compare data1 and data2 ...\n";
  double tol = 1.0e-20;
  bool r = CompareHierarchyData(d1, d2, tol);
  return r;
}
