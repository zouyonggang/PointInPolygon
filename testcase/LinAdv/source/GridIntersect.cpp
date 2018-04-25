#ifndef GRID_INTERSECT_CPP
#define GRID_INTERSECT_CPP

#include "GridIntersect.h"

#include <assert.h>
#include <string>

struct PointDescriptor {
  double x;
  double y;
#if (DIM == 3)
  double z;
#endif
};

struct CellItemDescriptor {
  std::vector<PointDescriptor> point;
}

template <int DIM>
class GridIntersectImpl {
public:
  GridIntersectImpl(tbox::Pointer<hier::Patch<DIM>> dest_patch) {
    int cell_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::CELL, 0);
    int node_number = dest_patch->getNumberOfEntities(
        hier::EntityUtilities::EntityType::NODE, 0);

    tbox::Pointer<pdat::NodeData<NDIM, double>> node =
        patch->getPatchGeometry->getNodeCoordinates();
    //确保坐标的深度为1，群组为维度
    assert(node->getDepth() == 1);
    assert(node->getGroup() == DIM);
    const double* node_coordinates = node->getPointer(0);

    //获取topology相关信息
    tbox::Pointer<hier::PatchTopology<NDIM>> topology =
        patch->getPatchTopology();
    tbox::Array<int> cell_adj_nodes_extent;
    tbox::Array<int> cell_adj_nodes_indices;
    topology->getCellAdjacencyNodes(cell_adj_nodes_extent,
                                    cell_adj_nodes_indices);
  }

  ~GridIntersectImpl() {}

  std::vector<int> PointInGrid(const double* points) {
    std::vector<int> temp;
    temp.push_back(1);
    return temp;
  }

  void RayIntersectGrid(const double* start_points, const double* direction,
                        std::vector<int> ids,
                        double* intersection_coordinates) {}

  void GridIntersectGrid(tbox::Pointer<hier::Patch<DIM>> src_patch, int number,
                         std::vector<int> ids) {}

private:
  std::vector<CellItemDescriptor> dest_patch_;  //比对目的网格片
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
std::vector<int> GridIntersect<DIM>::PointInGrid(const double* points) {
  return impl_->PointInGrid(points);
}

template <int DIM>
void GridIntersect<DIM>::RayIntersectGrid(const double* start_points,
                                          const double* direction,
                                          std::vector<int> ids,
                                          double* intersection_coordinates) {
  impl_->RayIntersectGrid(start_points, direction, ids,
                          intersection_coordinates);
}

template <int DIM>
void GridIntersect<DIM>::GridIntersectGrid(
    tbox::Pointer<hier::Patch<DIM>> src_patch, int number,
    std::vector<int> ids) {
  impl_->GridIntersectGrid(src_patch, number, ids);
}
#endif  // GRID_INTERSECT_CPP