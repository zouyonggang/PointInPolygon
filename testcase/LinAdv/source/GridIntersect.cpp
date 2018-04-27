#ifndef GRID_INTERSECT_CPP
#define GRID_INTERSECT_CPP

#include "GridIntersect.h"

#include "BoundingBox.h"
#include "DoubleVector.h"
#include "IntervalTree.h"
#include "NodeData.h"
#include "PatchGeometry.h"
#include "PatchTopology.h"

#include <assert.h>
#include <string>

struct BoxSpec {
  std::string id;
  int dim;
  std::string attr;

  BoxSpec(const std::string& id_, int dim_, const std::string& attr_)
      : id(id_), dim(dim_), attr(attr_) {}
};

template <int DIM>
class GridIntersectImpl {
public:
  struct CellItemDescriptor {
    int index;
    std::vector<hier::DoubleVector<DIM>> points;
    hier::BoundingBox<DIM> bbox;
  };

  GridIntersectImpl(tbox::Pointer<hier::Patch<DIM>> dest_patch) {
    int cell_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::CELL, 0);
    int node_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::NODE, 0);

    tbox::Pointer<pdat::NodeData<NDIM, double>> node =
        dest_patch->getPatchGeometry()->getNodeCoordinates();
    //确保坐标的深度为1，群组为维度
    assert(node->getDepth() == 1);
    assert(node->getGroup() == DIM);
    const double* nodes_coordinates = node->getPointer(0);

    //获取topology相关信息
    tbox::Pointer<hier::PatchTopology<NDIM>> topology =
        dest_patch->getPatchTopology();
    tbox::Array<int> cell_adj_nodes_extent;
    tbox::Array<int> cell_adj_nodes_indices;
    topology->getCellAdjacencyNodes(cell_adj_nodes_extent,
                                    cell_adj_nodes_indices);
    //获取每个网格单元节点坐标
    for (int i = 0; i < cell_number; i++) {  //遍历每个单元
      std::vector<hier::DoubleVector<DIM>> points;
      for (int node_iter = cell_adj_nodes_extent[i];
           node_iter < cell_adj_nodes_extent[i + 1];
           node_iter++) {  //遍历该单元的每个节点
        hier::DoubleVector<DIM> node_coordinate;
#if (DIM == 2)
        node_coordinate = hier::DoubleVector<DIM>(
            nodes_coordinates[cell_adj_nodes_indices[node_iter]],
            nodes_coordinates[cell_adj_nodes_indices[node_iter] + 1]);
#elif (DIM == 3)
        node_coordinate = hier::DoubleVector<DIM>(
            nodes_coordinates[cell_adj_nodes_indices[node_iter]],
            nodes_coordinates[cell_adj_nodes_indices[node_iter] + 1],
            nodes_coordinates[cell_adj_nodes_indices[node_iter] + 2]);
#endif
        points.push_back(node_coordinate);
      }

      //获取该网格单元的包围盒
      hier::BoundingBox<DIM> bbox;
      hier::DoubleVector<DIM> bbox_lower(1.0e16);
      hier::DoubleVector<DIM> bbox_upper(-1.0e16);
      for (int p = 0; p < points.size(); p++) {
        bbox_lower.min(points[p]);
        bbox_upper.max(points[p]);
      }
      bbox.setBoundingBox(bbox_lower, bbox_upper, i);
      //保存该网格单元
      CellItemDescriptor cell;
      cell.index = i;
      cell.points = points;
      cell.bbox = bbox;
      dest_patch_[i] = cell;
    }
  }

  ~GridIntersectImpl() {}

  /**
   * @brief
   * 输入一个点，将该点与patch内所有网格单元包围盒相交，输出与相交的网格单元索引
   *
   * @param points 输入参数，将要与网格单元相交的点坐标
   *
   * @return std::vector<int> 输出参数，与输入点相交的网格单元索引集合
   */
  std::vector<int> getBoundingBoxIntersect(const double* point) {
    //构造区间树，将patch的所有网格单元存入区间书的叶子节点
    hier::IntervalTree<DIM> interval_tree(dest_patch_.size());
    typename std::map<int, CellItemDescriptor>::iterator iter;
    for (; iter != dest_patch_.end(); iter++) {
      CellItemDescriptor& cell = iter->second;
      hier::DoubleVector<DIM> lower = cell.bbox.getLower();
      hier::DoubleVector<DIM> upper = cell.bbox.getUpper();
      double bbox[DIM * 2];
      for (int i = 0; i < DIM; i++) {
        bbox[i] = lower[i];
        bbox[DIM + i] = upper[i];
      }
      interval_tree.addElement(iter->first, bbox);
    }
    std::vector<int> result;
    interval_tree.getElementsListFromRange(point, point, result);
    return result;
  }

  std::vector<int> pointInGrid(const double* points, int n) {
    std::vector<int> result;
    for (int i = 0; i < n; i++) {
      double point[DIM];
      for (int d = 0; d < IDM; d++) point[d] = points[i * DIM + d];
      //第一步，通过包围盒和区间树大致判断点在那些单元之间，获取这些单元集合
      std::vector<int> bbox_intersect_cells = getBoundingBoxIntersect(point);
      if (bbox_intersect_cells.size() == 0) return;

      //第二步，遍历第一步得到的单元集合
      for (int c = 0; c < bbox_intersect_cells.size(); c++) {
        CellItemDescriptor cell = dest_patch_[bbox_intersect_cells[c]];

        //遍历该单元的所有face，从点出发取x方向的射线，利用rayIntersectGrid2计算是否与face相交，
        //根据相交情况确定crossing number，根据crossing
        // number的奇偶性来判断点是否在单元内
      }
    }
  }

  void rayIntersectGrid(const double* start_points, const double* direction,
                        std::vector<int> ids,
                        double* intersection_coordinates) {}

  void gridIntersectGrid(tbox::Pointer<hier::Patch<DIM>> src_patch, int number,
                         std::vector<int> ids) {}

private:
  std::map<int, CellItemDescriptor> dest_patch_;  //目的网格片
};

//
// 指向GridIntersect的实现
//
template <int DIM>
GridIntersect<DIM>::GridIntersect(tbox::Pointer<hier::Patch<DIM>> dest_patch)
    : impl_(new GridIntersectImpl<DIM>(dest_patch)) {}

template <int DIM>
GridIntersect<DIM>::~GridIntersect() {}

template <int DIM>
std::vector<int> GridIntersect<DIM>::pointInGrid(const double* points, int n) {
  return impl_->pointInGrid(points, n);
}

template <int DIM>
void GridIntersect<DIM>::rayIntersectGrid(const double* start_points,
                                          const double* direction, int n,
                                          std::vector<int> ids,
                                          double* intersection_coordinates) {
  impl_->rayIntersectGrid(start_points, direction, n, ids,
                          intersection_coordinates);
}

template <int DIM>
void GridIntersect<DIM>::gridIntersectGrid(
    tbox::Pointer<hier::Patch<DIM>> src_patch, int number,
    std::vector<int> ids) {
  impl_->gridIntersectGrid(src_patch, number, ids);
}
#endif  // GRID_INTERSECT_CPP