//
// 文件名: main.C
// 软件包: JAUMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 136 $
// 修改  : $Date: 2011-07-22 08:41:30 +0800 (五, 2011-07-22) $
// 描述  : 主控程序(矩形非结构网格上, 求解线性对流对流问题).
//

#include <string>
#include <vector>
using namespace std;

#include "Array.h"
#include "ArrayData.h"
#include "BoundingBox.h"
#include "CompareData.h"
#include "DoubleVector.h"
#include "EntityUtilities.h"
#include "GridGeometry.h"
#include "GridTopology.h"
#include "HierarchyTimeIntegrator.h"
#include "InputManager.h"
#include "IntervalTree.h"
#include "JAUMINManager.h"
#include "JAUMIN_config.h"
#include "JaVisDataWriter.h"
#include "PatchHierarchy.h"
#include "RestartManager.h"
#include "TimerManager.h"
#include "Utilities.h"
#include "VariableDatabase.h"
using namespace JAUMIN;

#include "LinAdv.h"
#include "LinAdvLevelIntegrator.h"

/*!
*************************************************************************
* @brief   在网格文件名前面追加input文件的路径
************************************************************************
*/
static void prefixInputDirName(const string& input_filename,
                               tbox::Pointer<tbox::Database> input_db) {
  string path_name = "";
  string::size_type slash_pos = input_filename.find_last_of('/');
  if (slash_pos != string::npos)
    path_name = input_filename.substr(0, slash_pos + 1);

  string mesh_file = input_db->getDatabase("GridGeometry")
                         ->getDatabase("MeshImportationParameter")
                         ->getString("file_name");

  slash_pos = mesh_file.find_first_of('/');
  if (slash_pos != 0) {
    input_db->getDatabase("GridGeometry")
        ->getDatabase("MeshImportationParameter")
        ->putString("file_name", path_name + mesh_file);
  }
}

struct PointDescriptor {
  std::vector<double> data;
};

/*!
*************************************************************************
*
* @brief 基于JAUMIN框架的非结构网格, 求解线性对流问题.
*
* 该程序分以下几个步骤:
* -# 预处理: 初始化MPI和JAUMIN环境, 解析输入文件, 读取主程序控制参数;
* -# 创建网格片层次结构和时间积分算法类对象, 主要包括:
*       - 网格几何       hier::GridGeometry<NDIM>
*       - 网格拓扑       hier::GridTopology<NDIM>
*       - 网格片层次结构 hier::PatchHierarchy<NDIM>
*    -# 网格片积分算法 LinAdv
*       - 应用级: 提供求解线性对流方程的数值计算子程序
*    -# 网格层积分算法 LinAdvLevelIntegrator
*       - 应用级: 提供基于网格层的线性对流方程求解流程
*    -# 网格层时间积分算法 algs::HierarchyTimeIntegrator<NDIM>
* -# 初始化网格片层次结构和物理量数据片;
* -# 循环: 时间步进;
* -# 后处理: 释放应用类对象, 释放JAUMIN和MPI内部资源.
*
************************************************************************
*/
int main(int argc, char* argv[]) {
  // 初始化MPI和JAUMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JAUMINManager::startup();

  {
    /*******************************************************************************
     *                               预  处  理 *
     *******************************************************************************/
    // 解析命令行参数:
    string input_filename;  // 输入文件名.
    if (argc != 2) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                 << "<restart dir> <restore number> [options]\n"
                 << "  options:\n"
                 << "  none at this time" << endl;
      tbox::MPI::abort();
      return (-1);
    } else {
      input_filename = argv[1];
    }

    /// 把信息输出到log文件
    tbox::plog << "input_filename = " << input_filename << endl;

    // 解析输入文件的计算参数到输入数据库, 称之为根数据库.
    tbox::Pointer<tbox::Database> input_db =
        new tbox::InputDatabase("input_db");
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    // 在网格文件名前面追加input文件的路径
    prefixInputDirName(input_filename, input_db);

    // 从根数据库中获得名称为"Main"的子数据库.
    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象 *
     *******************************************************************************/
    //(1) 创建非结构网格几何.
    tbox::Pointer<hier::GridGeometry<NDIM>> grid_geometry =
        new hier::GridGeometry<NDIM>("GridGeometry",
                                     input_db->getDatabase("GridGeometry"));

    //(2) 创建非结构网格拓扑.
    tbox::Pointer<hier::GridTopology<NDIM>> grid_topology =
        new hier::GridTopology<NDIM>("GridTopology",
                                     input_db->getDatabase("GridTopology"));

    //(3) 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry,
                                       grid_topology, true);

    //(4) 创建网格片时间积分算法类（应用级:
    //提供求解线性对流方程的数值计算子程序）.
    LinAdv* linadv_advection_model =
        new LinAdv("LinAdv", input_db->getDatabase("LinAdv"));

    //(5) 创建网格层时间积分算法类（应用级:
    //提供基于网格层的线性对流方程求解流程）.
    LinAdvLevelIntegrator* linadv_level_integrator = new LinAdvLevelIntegrator(
        "LinAdvLevelIntegrator", linadv_advection_model);

    //(6) 创建网格层时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM>> time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            linadv_level_integrator, true);

    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    /************************************************************************************
     *                              输 出 数 据 *
     ************************************************************************************/
    {
      //获取网格层
      int number_levels = patch_hierarchy->getNumberOfLevels();
      cout << endl << "网格层总数：" << number_levels << endl;
      for (int level_id = 0; level_id < number_levels; level_id++) {
        cout << "当前网格层：" << level_id << endl;
        tbox::Pointer<hier::PatchLevel<NDIM>> patch_level =
            patch_hierarchy->getPatchLevel(level_id);

        for (typename hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++) {
          tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++"
                     << endl;
          tbox::Pointer<hier::Patch<NDIM>> patch = patch_level->getPatch(p());
          cout << "patch index:" << patch->getIndex() << endl;
          int cell_number = patch->getNumberOfEntities(
              hier::EntityUtilities::EntityType::CELL, 2);
          cout << "Patch cell entity:" << cell_number << endl;
          int node_number = patch->getNumberOfNodes(2);
          cout << "Patch node entity:" << node_number << endl << endl;

          //获取geometry相关信息
          tbox::Pointer<hier::PatchGeometry<NDIM>> geometry =
              patch->getPatchGeometry();
          cout << "the number of geometry entity set:"
               << geometry->getNumberOfEntitySet() << endl;
          tbox::Pointer<pdat::CellData<NDIM, double>> cell_coordinates =
              geometry->getCellCoordinates();
          cout << "cell coordinates depth：" << cell_coordinates->getDepth()
               << endl;
          cout << "cell coordinates groups：" << cell_coordinates->getGroup()
               << endl;
          double* cell_coordinate = cell_coordinates->getPointer();
          cout << "cell coordinates：";
          for (int i = 0; i < cell_number * 2; i++) {
            cout << cell_coordinate[i] << " ";
          }
          tbox::Pointer<pdat::NodeData<NDIM, double>> node_coordinates =
              geometry->getNodeCoordinates();
          cout << "node coordinates depth：" << node_coordinates->getDepth()
               << endl;
          cout << "node coordinates groups：" << node_coordinates->getGroup()
               << endl;
          const double* coordinates = node_coordinates->getPointer(0);
          cout << "node coordinates：";
          for (int i = 0; i < node_number; i++) {
            for (int j = 0; j < NDIM; j++)
              cout << coordinates[i * NDIM + j] << " ";
            cout << ",";
          }
          cout << endl;

          //获取bounding box
          cout << "bounding box:";
          // for (int i = 0; i <= node_number; i++) {
          // cout << "  bounding box entity " << i << ":";
          // tbox::Array<int> node_flag(4, 10);
          // hier::BoundingBox<NDIM> bbox = geometry->getBoundingBox(node_flag);
          hier::DoubleVector<NDIM> a(0, 0);
          hier::DoubleVector<NDIM> b(5, 2);
          hier::BoundingBox<NDIM> bbox(a, b, 10);
          cout << "index" << bbox.getIndex() << endl;
          hier::DoubleVector<NDIM> lower = bbox.getLower();
          hier::DoubleVector<NDIM> upper = bbox.getUpper();
          cout << "bbox coordiante lower:" << lower(0) << "," << lower(1)
               << " ";
          cout << "bbox coordiante upper:" << upper[0] << "," << upper[1]
               << endl;

          //区间树
          double box1[4] = {0, 5, 0, 2};
          double box2[4] = {7, 10, 0, 2};
          double test_box_lo[4] = {7, 1};
          double test_box_up[4] = {8, 2};

          hier::IntervalTree<NDIM> interval_tree(2);
          interval_tree.addElement(1, box1);
          interval_tree.addElement(2, box2);
          interval_tree.constructTree();

          std::vector<int> result;
          interval_tree.getElementsListFromRange(test_box_lo, test_box_up,
                                                 result);
          cout << "result size:" << result.size() << endl;
          for (int i = 0; i < result.size(); i++) cout << result[i] << " ";

        }  // end patch
      }    // end level
    }      // end cout
  }

  // 释放JAUMIN和MPI内部资源.
  tbox::JAUMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
