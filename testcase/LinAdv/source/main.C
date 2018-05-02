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

#include "ArrayData.h"
#include "CompareData.h"
#include "EntityUtilities.h"
#include "GridGeometry.h"
#include "GridTopology.h"
#include "HierarchyTimeIntegrator.h"
#include "InputManager.h"
#include "JAUMINManager.h"
#include "JAUMIN_config.h"
#include "JaVisDataWriter.h"
#include "PatchHierarchy.h"
#include "RestartManager.h"
#include "TimerManager.h"
#include "Utilities.h"
#include "VariableDatabase.h"
using namespace JAUMIN;

#include "GridIntersect.h"
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
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geometry =
        new hier::GridGeometry<NDIM>("GridGeometry",
                                     input_db->getDatabase("GridGeometry"));

    //(2) 创建非结构网格拓扑.
    tbox::Pointer<hier::GridTopology<NDIM> > grid_topology =
        new hier::GridTopology<NDIM>("GridTopology",
                                     input_db->getDatabase("GridTopology"));

    //(3) 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
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
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            linadv_level_integrator, true);

    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    /************************************************************************************
     *                              测 试 数 据 *
     ************************************************************************************/
    {
      //获取网格层
      int number_levels = patch_hierarchy->getNumberOfLevels();
      cout << endl << "网格层总数：" << number_levels << endl;
      for (int level_id = 0; level_id < number_levels; level_id++) {
        cout << "当前网格层：" << level_id << endl;
        tbox::Pointer<hier::PatchLevel<NDIM> > patch_level =
            patch_hierarchy->getPatchLevel(level_id);

        for (typename hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++) {
          //存储网格内部点
          std::vector<PointDescriptor> inside_points;
          //存储边界点
          std::vector<PointDescriptor> boundary_points;
          //存储网格外部点
          std::vector<PointDescriptor> outside_points;
          tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++"
                     << endl;
          tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(p());
          cout << "patch index:" << patch->getIndex() << endl;
          int cell_number =
              patch->getNumberOfEntities(hier::EntityUtilities::CELL, 0);
          cout << "Patch cell entity:" << cell_number << endl;
          int node_number =
              patch->getNumberOfEntities(hier::EntityUtilities::NODE, 0);
          tbox::Pointer<hier::PatchGeometry<NDIM> > geometry =
              patch->getPatchGeometry();
          tbox::Pointer<pdat::CellData<NDIM, double> > cell_coordinates =
              geometry->getCellCoordinates();
          double* cell_coordinate = cell_coordinates->getPointer();
          tbox::Pointer<pdat::NodeData<NDIM, double> > node_coordinates =
              geometry->getNodeCoordinates();
          const double* node_coordinate = node_coordinates->getPointer(0);

          //将网格中心点赋值给inside-points
          for (int i = 0; i < cell_number; i++) {
            PointDescriptor point;
            for (int j = 0; j < NDIM; j++)
              point.data.push_back(cell_coordinate[j + i * NDIM]);
            inside_points.push_back(point);
          }
          //将x轴最下方的顶点赋值给boundary-points
          PointDescriptor x_lower_point = inside_points[0];
          for (int i = 0; i < node_number; i++)
            if (node_coordinate[i * NDIM] < x_lower_point.data[0])
              for (int j = 0; j < NDIM; j++)
                x_lower_point.data[j] = node_coordinate[i * NDIM + j];
          boundary_points.push_back(x_lower_point);
          //将内部影像区顶点赋值给outside_points
          int node_number_with_ghost = patch->getNumberOfNodes(2);
          for (int i = node_number; i < node_number_with_ghost; i++) {
            PointDescriptor point;
            for (int j = 0; j < NDIM; j++)
              point.data.push_back(node_coordinate[i * NDIM + j]);
            outside_points.push_back(point);
          }

          //输出inside_points,boundary_points,outside_points
          cout << "inside_point:";
          double inside_point[inside_points.size() * NDIM];
          int iter = 0;
          for (int i = 0; i < inside_points.size(); i++) {
            for (int j = 0; j < NDIM; j++) {
              inside_point[iter++] = inside_points[i].data[j];
              cout << inside_points[i].data[j] << " ";
            }
            cout << ";";
          }

          cout << endl << "boundary_points:";
          double boundary_point[boundary_points.size() * NDIM];
          iter = 0;
          for (int i = 0; i < boundary_points.size(); i++) {
            for (int j = 0; j < NDIM; j++) {
              boundary_point[iter++] = boundary_points[i].data[j];
              cout << boundary_points[i].data[j] << " ";
            }
            cout << ";";
          }
          cout << endl << "outside_points:";
          double outside_point[outside_points.size() * NDIM];
          iter = 0;
          for (int i = 0; i < outside_points.size(); i++) {
            for (int j = 0; j < NDIM; j++) {
              outside_point[iter++] = outside_points[i].data[j];
              cout << outside_points[i].data[j] << " ";
            }
            cout << ";";
          }
          cout << endl;

          //进行点与网格定位
          tbox::Pointer<GridIntersect> intersect = new GridIntersect(patch);
          vector<int> inside_result =
              intersect->pointInGrid(inside_point, inside_points.size());
          cout << "inside相交网格编号为:";
          for (int i = 0; i < inside_result.size(); i++)
            cout << inside_result[i] << " ";

          vector<int> outside_result =
              intersect->pointInGrid(outside_point, outside_points.size());
          cout << endl << "outside相交网格编号为:";
          for (int i = 0; i < outside_result.size(); i++)
            cout << outside_result[i] << " ";

          // vector<int> boundary_result =
          //     intersect->pointInGrid(boundary_point, boundary_points.size());
          // cout << endl << "boundary相交网格编号为:";
          // for (int i = 0; i < boundary_result.size(); i++)
          //   cout << boundary_result[i] << " ";

        }  // end patch
      }    // end level
    }      // end test
  }

  // 释放JAUMIN和MPI内部资源.
  tbox::JAUMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
