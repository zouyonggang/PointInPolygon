#ifndef GRID_INTERSECT_CPP
#define GRID_INTERSECT_CPP

#include "GridIntersect.h"

#include "Array.h"
#include "BoundingBox.h"
#include "DoubleVector.h"
#include "IntervalTree.h"
#include "NodeData.h"
#include "PatchGeometry.h"
#include "PatchTopology.h"

#include <assert.h>
#include <string>

/**
 * @brief 描述一个3维结点，里面包含patch内节点编号，节点坐标
 *
 */
struct NodeDescriptor {
  int index;
  hier::DoubleVector<3> coord;
};

/**
 * @brief 描述一个3维棱，包含patch内棱编号，棱周围节点编号
 *
 */
struct EdgeDescriptor {
  int index;
  std::vector<int> nodes;
};

/**
 * @brief 描述一个面，包含patch内面编号，面周围节点编号、棱编号
 *
 */
struct FaceDescriptor {
  int index;
  std::vector<int> nodes;
  std::vector<int> edges;
  std::vector<int> cells;  //用来获取网格外表面
};

/**
 * @brief 描述一个三维体，包含patch体编号，体周围节点编号、棱编号、面编号
 *
 */
struct CellDescriptor {
  int index;
  std::vector<int> nodes;
  std::vector<int> edges;
  std::vector<int> faces;
};

class GridIntersectImpl {
public:
  /**
   * @brief构造函数，从jaumin框架内得到点、边、面、体（非影像区），并保存到该类结构中
   *
   * @param dest_patch 输入的patch网格片
   */
  GridIntersectImpl(tbox::Pointer<hier::Patch<3>> dest_patch) {
    int node_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::NODE, 0);
    int edge_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::EDGE, 0);
    int face_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::FACE, 0);
    int cell_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::CELL, 0);
    tbox::Pointer<hier::PatchTopology<3>> topology =
        dest_patch->getPatchTopology();

    //获取节点,保存节点
    tbox::Pointer<hier::PatchGeometry<3>> patch_geometry =
        dest_patch->getPatchGeometry();
    tbox::Pointer<pdat::NodeData<3, double>> node =
        patch_geometry->getNodeCoordinates();
    assert(node->getDepth() == 1);
    assert(node->getGroup() == 3);
    const double* nodes_coordinates = node->getPointer(0);
    for (int i = 0; i < node_number; i++) {
      NodeDescriptor node;
      node.index = i;
      node.coord = hier::DoubleVector<3>(nodes_coordinates[i * 3],
                                         nodes_coordinates[i * 3 + 1],
                                         nodes_coordinates[i * 3 + 2]);
      patch_nodes_[i] = node;
    }
    //获取边
    tbox::Array<int> edge_adj_nodes_extent;
    tbox::Array<int> edge_adj_nodes_indices;
    topology->getEdgeAdjacencyNodes(edge_adj_nodes_extent,
                                    edge_adj_nodes_indices);
    for (int i = 0; i < edge_number; i++) {
      EdgeDescriptor edge;
      edge.index = i;
      for (int j = edge_adj_nodes_extent[i]; j < edge_adj_nodes_extent[i + 1];
           j++)
        edge.nodes.push_back(edge_adj_nodes_indices[j]);
      patch_edges_[i] = edge;
    }

    //获取面
    tbox::Array<int> face_adj_nodes_extent;
    tbox::Array<int> face_adj_nodes_indices;
    topology->getFaceAdjacencyNodes(face_adj_nodes_extent,
                                    face_adj_nodes_indices);
    tbox::Array<int> face_adj_edges_extent;
    tbox::Array<int> face_adj_edges_indices;
    topology->getFaceAdjacencyEdges(face_adj_edges_extent,
                                    face_adj_edges_indices);
    tbox::Array<int> face_adj_cells_extent;
    tbox::Array<int> face_adj_cells_indices;
    topology->getFaceAdjacencyCells(face_adj_cells_extent,
                                    face_adj_cells_indices);
    for (int i = 0; i < face_number; i++) {
      FaceDescriptor face;
      face.index = i;
      for (int j = face_adj_nodes_extent[i]; j < face_adj_nodes_extent[i + 1];
           j++)
        face.nodes.push_back(face_adj_nodes_indices[j]);
      for (int j = face_adj_edges_extent[i]; j < face_adj_edges_extent[i + 1];
           j++)
        face.edges.push_back(face_adj_edges_indices[j]);
      //获取面周围的体，并判断该面是否为外表面
      for (int j = face_adj_cells_extent[i]; j < face_adj_cells_extent[i + 1];
           j++) {
        if (face_adj_cells_extent[i + 1] - face_adj_cells_extent[i] == 1)
          patch_outer_faces_[i] = face_adj_cells_indices[j];
        face.cells.push_back(face_adj_cells_indices[j]);
      }
      patch_faces_[i] = face;
    }

    //获取体
    tbox::Array<int> cell_adj_nodes_extent;
    tbox::Array<int> cell_adj_nodes_indices;
    topology->getCellAdjacencyNodes(cell_adj_nodes_extent,
                                    cell_adj_nodes_indices);
    tbox::Array<int> cell_adj_edges_extent;
    tbox::Array<int> cell_adj_edges_indices;
    topology->getCellAdjacencyEdges(cell_adj_edges_extent,
                                    cell_adj_edges_indices);
    tbox::Array<int> cell_adj_faces_extent;
    tbox::Array<int> cell_adj_faces_indices;
    topology->getCellAdjacencyFaces(cell_adj_faces_extent,
                                    cell_adj_faces_indices);
    for (int i = 0; i < cell_number; i++) {
      CellDescriptor cell;
      cell.index = i;
      for (int j = cell_adj_nodes_extent[i]; j < cell_adj_nodes_extent[i + 1];
           j++)
        cell.nodes.push_back(cell_adj_nodes_extent[j]);
      for (int j = cell_adj_edges_extent[i]; j < cell_adj_edges_extent[i + 1];
           j++)
        cell.edges.push_back(face_adj_edges_indices[j]);
      for (int j = cell_adj_faces_extent[i]; j < cell_adj_faces_extent[i + 1];
           j++)
        cell.faces.push_back(cell_adj_faces_indices[j]);
      patch_cells_[i] = cell;
    }
  }

  ~GridIntersectImpl() {}

  /**
   * @brief 判断点2，是否在点0，点1组成的线段的左边，只能适用于2维情况
   *
   * @param point0 点0，长度为2
   * @param point1 点1，长度为2
   * @param point2 点2，长度为2
   * @return int >0 p2 在通过p0,p1线段的左边
   *             =0 p2 在通过p0,p1线段上
   *             <0 p2 在通过p0,p1线段的右边
   */
  inline int isLeft(const double* point0, const double* point1,
                    const double* point2) {
    assert(point1 != NULL && point2 != NULL && point0 != NULL);
    return ((point1[0] - point0[0]) * (point2[1] - point0[1]) -
            (point2[0] - point0[0]) * (point1[1] - point0[1]));
  }

  /**
   * @brief 利用缠绕算法，计算二维平面内，点是否在多边形内
   *
   * @param point 输入点
   * @param vertex 输入多边形的顶点
   * @param n 多边形顶点的个数
   * @return int >0 点在多边形内部
   *             =0 点在多边形外部
   *             <0 算法错误情况
   */
  int wnPnpoly(const double* point, const double* vertex, int n) {
    assert(point != NULL && vertex != NULL);
    assert(n > 1);

    int wn = 0;  //记录winding number 数量

    //遍历多边形的所有边
    for (int i = 0; i < n; i++) {  // vertex[i]和vertex[i+1]组成的边
      double V[2] = {vertex[i * 2], vertex[i * 2 + 1]};
      double V_next[2] = {vertex[(i + 1) * 2], vertex[(i + 1) * 2 + 1]};
      if (V[1] <= point[1]) {                //输入点在线段下顶点上方
        if (V_next[1] > point[1])            //一个向上的缠绕
          if (isLeft(V, V_next, point) > 0)  // p在边的左边
            ++wn;

      } else {                               //输入点在线段下顶点上方
        if (V_next[1] <= point[1])           //一个向下的缠绕
          if (isLeft(V, V_next, point) < 0)  // p在边的右边
            --wn;
      }
    }
    assert(wn >= 0);
    return wn;
  }

  /**
   * @brief 计算射线与平面的交点
   *
   * @param point 射线起点
   * @param direction 射线方向
   * @param face 平面（由三个点构成）
   * @param intersection 输出参数，交点坐标
   * @param dist 输出参数，射线起点到交点距离
   * @param tag 输出参数，出错信息，
   *              0：成功求出交点
   *             -1：射线与平面平行
   *             -2：射线所在直线与平面的交点在射线的副方向
   */
  void lineSurfaceIntersect(const double* point, const double* direction,
                            const double* face, double* intersection,
                            double& dist, int tag) {
    assert(point != NULL && direction != NULL && face != NULL);

    const double epson = 1.0e-12;
    dist = -1.0e-15;
    double x1, y1, z1, x2, y2, z2, AA, BB, CC, DD;
    double delta, tmp, t;

    //平面上的两个切方向
    x1 = face[1] - face[0];
    y1 = face[1 * 3 + 1] - face[1 * 3];
    z1 = face[2 * 3 + 1] - face[2 * 3];
    x2 = face[2] - face[0];
    y2 = face[1 * 3 + 2] - face[1 * 3];
    z2 = face[2 * 3 + 2] - face[2 * 3];
    //平面法向，平面方程 AA x + BB y + CC z + DD = 0
    AA = y1 * z2 - z1 * y2;
    BB = z1 * x2 - x1 * z2;
    CC = x1 * y2 - y1 * x2;
    DD = -1 * face[0] * AA - face[1 * 3] * BB - face[2 * 3] * CC;
    //平面发向与射线切向的内积
    delta = AA * direction[0] + BB * direction[1] + CC * direction[2];
    //射线方程(xm代表射线起点，xd代表射线方向) x = xm + xd t
    //                                   y = ym + yd t
    //                                   z = zm + zd t
    if (abs(delta) <= epson)
      //射线切向与平面平行
      tag = -1;
    else {
      //射线与平面不平行，计算交点对应的t值
      tmp = -1 * AA * point[0] - BB * point[1] - CC * point[2] - DD;
      t = tmp / delta;
      double compare = 0;
      if (t >= compare)  // t>0,解有意义，计算交点坐标
      {
        tag = 0;
        intersection[0] = point[0] + direction[0] * t;
        intersection[1] = point[1] + direction[1] * t;
        intersection[2] = point[2] + direction[2] * t;
        dist = t * sqrt(intersection[0] * intersection[0] +
                        intersection[1] * intersection[1] +
                        intersection[2] * intersection[2]);
      } else
        // t < 0，说明射线与平面的交点在射线的负方向，射线本上与平面五交点
        tag = -2;
    }
  }

  /**
   * @brief 将同一个面三维点降为二维点
   *
   * @param points 输入三维点集合
   * @param n 点个数
   * @param result 返回二维点集合
   */
  void reduceDim(double* points, int n, double* result) {
    int dim_index = 0;
    for (int d = 0; d < 3; d++) {
      int i = 0;
      for (; i < n; i++)
        if (points[i * 3 + d] != points[(i + 1) * 3 + d]) break;
      assert(i == n);  // points not in a face
      dim_index = d;
    }
    for (int i = 0; i < n; i++) {
      if (dim_index == 0) {
        result[i * 2] = points[i * 3 + 1];
        result[i * 2 + 1] = points[i * 3 + 2];
      } else if (dim_index == 1) {
        result[i * 2] = points[i * 3];
        result[i * 2 + 1] = points[i * 3 + 2];
      } else {
        result[i * 2] = points[i * 3];
        result[i * 2 + 1] = points[i * 3 + 1];
      }
    }
  }

  /**
   * @brief 计算射线与点组成的面相交，并判断交点是否在组成面的多边形中
   *
   * @param point 射线起点
   * @param direction 射线方向
   * @param face 组成面的点
   * @param n 组成面的点的个数
   * @param intersection 输出参数，交点坐标
   * @param dist 输出参数，射线起点到交点的距离
   * @param tag 输出参数，标识符
   *              0：成功求出交点，并且交点在face面内部
   *             -1：射线与平面平行
   *             -2：射线所在直线与平面的交点在射线的副方向
   *             -3：射线与所在平面相交，但是交点在face外
   */
  void rayIntersectFace(const double* point, const double* direction,
                        const double* face, int n, double* intersection,
                        double& dist, int tag) {
    assert(n >= 3);
    //取前三个点组成的面与射线相交，获取交点与距离
    double sub_face[9];
    for (int i = 0; i < 9; i++) sub_face[i] = face[i];
    lineSurfaceIntersect(point, direction, sub_face, intersection, dist, tag);

    if (tag == 0) {  //有交点
      //将交点于组成面的点降为2维
      double new_points[(n + 1) * 3];
      for (int i = 0; i < n * 3; i++) new_points[i] = face[i];
      new_points[n * 3] = intersection[0];
      new_points[n * 3 + 1] = intersection[1];
      new_points[n * 3 + 2] = intersection[2];
      double result_points[(n + 1) * 2];
      reduceDim(new_points, n + 1, result_points);

      //判断交点是否在face构成的多边形中
      double new_intersection[2] = {result_points[n * 2],
                                    result_points[n * 2 + 1]};
      if (wnPnpoly(new_intersection, result_points, n) > 0)
        return;
      else
        tag == -3;
    }
  }

  /**
   * @brief
   * 输入一个点，将该点与patch内所有网格单元包围盒相交，输出与相交的网格单元索引
   *
   * @param points 输入参数，将要与网格单元相交的点坐标
   *
   * @return std::vector<int> 输出参数，与输入点相交的网格单元索引集合
   */
  std::vector<int> getBoundingBoxIntersect(const double* point) {
    assert(point != NULL);

    //构造区间树，将patch的所有网格单元bbox存入区间书的叶子节点
    hier::IntervalTree<3> interval_tree(patch_cells_.size());
    for (std::map<int, CellDescriptor>::iterator it = patch_cells_.begin();
         it != patch_cells_.end(); it++) {
      CellDescriptor& cell = it->second;
      double bbox[6] = {1.0e16, 1.0e16, 1.0e16, -1.0e16, -1.0e16, -1.0e16};
      for (int i = 0; i < cell.nodes.size(); i++) {
        assert(cell.nodes[i] <= patch_nodes_.size());
        NodeDescriptor& node = patch_nodes_[cell.nodes[i]];
        for (int d = 0; d < 3; d++) {
          bbox[d] = std::min(bbox[d], node.coord[d]);
          bbox[d + 3] = std::max(bbox[d + 3], node.coord[d]);
        }
      }
      interval_tree.addElement(it->first, bbox);
    }
    std::vector<int> result;
    interval_tree.getElementsListFromRange(point, point, result);
    return result;
  }

  /**
   * @brief 判断点是否在网格中
   *
   * @param points 输入点坐标
   * @param n 点的数量
   * @return std::vector<int> 返回每个点所在网格单元，
   *                          -1代表在网格外，
   *                          >=0代表在网格中
   */
  std::vector<int> pointInGrid(const double* points, int n) {
    assert(points != NULL);
    assert(n > 0);

    std::vector<int> result;
    for (int i = 0; i < n; i++) {
      double point[3];
      for (int d = 0; d < 3; d++) point[d] = points[i * 3 + d];
      //第一步，通过包围盒和区间树大致判断点在那些单元之间，获取这些单元集合
      std::vector<int> bbox_intersect_cells = getBoundingBoxIntersect(point);
      if (bbox_intersect_cells.size() == 0) return result;

      //第二步，遍历第一步得到的单元集合
      for (int c = 0; c < bbox_intersect_cells.size(); c++) {
        CellDescriptor cell = patch_cells_[bbox_intersect_cells[c]];

        //采用cross number的方法来判断点是否在网格单元中有一个弊端，
        //就是当点在单元边缘时,由于double误差关系，会使计算结果不准确，
        //所以取点上下左右前后附近六个点进行测试，如果这六个点都在网格内，
        //说明该点在网格内，如果有的在网格内，有的在网格外进行报错
        bool result[6] = {false};
        double around_points[18] = {
            point[0] + 1.0e-8, point[1],          point[2],
            point[0] - 1.0e-8, point[1],          point[2],
            point[0],          point[1] + 1.0e-8, point[2],
            point[0],          point[1] - 1.0e-8, point[2],
            point[0],          point[1],          point[2] + 1.0e-8,
            point[0],          point[1],          point[2] - 1.0e-8};
        for (int p = 0; p < 6; p++) {
          double around_point[3] = {around_points[p * 3],
                                    around_points[p * 3 + 1],
                                    around_points[p * 3 + 2]};

          //第三步，遍历该单元的所有face，从点出发取x方向的射线，利用rayIntersectFace计算是否与face相交
          //记录交点个数，如果个数为偶数则该点在网格单元内，如果为奇数，则在网格单元外
          int cross_number = 0;
          for (int f = 0; f < cell.faces.size(); f++) {
            FaceDescriptor face = patch_faces_[cell.faces[f]];
            std::vector<int> face_nodes = face.nodes;
            double face_nodes_coord[face_nodes.size() * 3];
            double direction[3] = {1, 0, 0};
            double intersection[3];
            double dist;
            int tag;
            rayIntersectFace(around_point, direction, face_nodes_coord,
                             face_nodes.size(), intersection, dist, tag);
            if (tag == 0) cross_number++;
          }  // end loop cell face
          if (cross > 0 && (cross_number % 2 == 0)) reuslt[i] = true;
        }

        for (int p = 0; p < 6; p++)
          if (result[p] != result[p + 1]) {
          }

      }  // end loop cell
    }    // end loop points
  }

  /**
   * @brief 计算输入射线集合与网格的交点，输出交点坐标以及交点所在网格单元
   *
   * @param points 射线起点
   * @param directions 射线方向
   * @param n 射线数量
   * @param ids 每条射线与网格相交，交点所在网格单元索引，
   * 如果射线不相交,索引值为-1
   * @param intersections_coordinates 每条射线与网格相交，交点坐标，
   * 如果射线不相交，交点坐标为{-1.0e16, -1.0e16, -1.0e16}
   */
  void rayIntersectGrid(const double* points, const double* directions, int n,
                        std::vector<int> ids,
                        double* intersections_coordinates) {
    assert(points != NULL && directions != NULL &&
           intersections_coordinates != NULL);
    //遍历每一条射线
    for (int r = 0; r < n; r++) {
      double point[3] = {points[r * 3], points[r * 3 + 1], points[r * 3 + 2]};
      double direction[3] = {directions[r * 3], directions[r * 3 + 1],
                             directions[r * 3 + 2]};
      double intersection[3] = {-1.0e16, -1.0e16, -1.0e16};
      double dist = 1.0e16;
      int intersect_face_index = -1;

      //遍历每一个外表面
      for (std::map<int, int>::iterator iter = patch_outer_faces_.begin();
           iter != patch_outer_faces_.end(); iter++) {
        FaceDescriptor face = patch_faces_[iter->first];
        double sub_intersection[3];
        double sub_dist;
        int tag;
        std::vector<int> face_nodes = face.nodes;
        double face_nodes_coord[face_nodes.size() * 3];
        for (int i = 0; i < face_nodes.size(); i++) {
          assert(face_nodes[i] <= patch_nodes_.size());
          face_nodes_coord[i * 3] = patch_nodes_[face_nodes[i]].coord[0];
          face_nodes_coord[i * 3 + 1] = patch_nodes_[face_nodes[i]].coord[1];
          face_nodes_coord[i * 3 + 2] = patch_nodes_[face_nodes[i]].coord[2];
        }
        //获取射线与face面的交点
        rayIntersectFace(point, direction, face_nodes_coord, face_nodes.size(),
                         sub_intersection, sub_dist, tag);
        if (tag == 0 && sub_dist < dist) {
          for (int i = 0; i < 3; i++) intersection[i] = sub_intersection[i];
          dist = sub_dist;
          intersect_face_index = iter->first;
        }
      }  // end loop patch outer face
      ids.push_back(intersect_face_index);
      for (int i = 0; i < 3; i++)
        intersections_coordinates[r * 3 + i] = intersection[i];
    }  // end loop ray
  }

  void gridIntersectGrid(tbox::Pointer<hier::Patch<3>> src_patch, int number,
                         std::vector<int> ids) {}

private:
  std::map<int, NodeDescriptor> patch_nodes_;  // patch内节点集合
  std::map<int, EdgeDescriptor> patch_edges_;  // patch内棱集合
  std::map<int, FaceDescriptor> patch_faces_;  // patch内面集合
  std::map<int, CellDescriptor> patch_cells_;  // patch内体集合
  std::map<int, int> patch_outer_faces_;  // patch外表面集合<face-id,cell-id>
};

//
// 指向GridIntersect的实现
//

GridIntersect::GridIntersect(tbox::Pointer<hier::Patch<3>> dest_patch)
    : impl_(new GridIntersectImpl(dest_patch)) {}

GridIntersect::~GridIntersect() {}

std::vector<int> GridIntersect::pointInGrid(const double* points, int n) {
  return impl_->pointInGrid(points, n);
}

void GridIntersect::rayIntersectGrid(const double* start_points,
                                     const double* direction, int n,
                                     std::vector<int> ids,
                                     double* intersection_coordinates) {
  impl_->rayIntersectGrid(start_points, direction, n, ids,
                          intersection_coordinates);
}

void GridIntersect::gridIntersectGrid(tbox::Pointer<hier::Patch<3>> src_patch,
                                      int number, std::vector<int> ids) {
  impl_->gridIntersectGrid(src_patch, number, ids);
}
#endif  // GRID_INTERSECT_CPP