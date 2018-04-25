#ifndef GRID_INTERSECT_H
#define GRID_INTERSECT_H

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include "Patch.h"
using namespace JAUMIN;

template <int DIM>
class GridIntersectImpl;

/**
 * @brief  GridIntersect
 * 该模块实现三个功能：1.查询输入点是否在输入网格内，如果在输出点所在网格单元
 * 2.给定网格外射线与网格，输出射线与网格外表面相加网格单元以及交点坐标
 * 3.给定源数据网格和目标单元网格，输出两个网格相交的源单元个数，以及源单元在网格片中的索引号
 *
 */
template <int DIM>
class GridIntersect : public boost::noncopyable {
public:
  /**
   * @brief 构造函数
   *
   * @param dest_patch 用于点相交、射线相交、网格相交的目的网格
   */
  GridIntersect(tbox::Pointer<hier::Patch<DIM>> dest_patch);
  ~GridIntersect();

  /**
   * @brief 查询点是否在网格内，如果在输出所在网格单元编号
   *
   * @param points 将要查询的点集合
   * @return vector<int> 每个点所在的网格单元编号
   */
  std::vector<int> PointInGrid(const double* points);

  /**
   * @brief 点与网格外表面相交,输出相交网格单元编号与交点坐标
   *
   * @param start_points 射线起点集合
   * @param direction 射线方向集合
   * @param ids 输出参数，交点所在网格单元编号
   * @param intersection 输出参数，交点坐标
   */
  void RayIntersectGrid(const double* start_points, const double* direction,
                        std::vector<int> ids, double* intersection_coordinates);

  /**
   * @brief 两个网格相交，输出相交网格单元个数以及单元在源网格片中的索引
   *
   * @param src_patch 源网格
   * @param number 相交单元个数
   * @param ids 源网格中相交网格单元索引
   */
  void GridIntersectGrid(tbox::Pointer<hier::Patch<DIM>> src_patch, int number,
                         std::vector<int> ids);

private:
  boost::shared_ptr<GridIntersectImpl<DIM>> impl_;
};

#include "GridIntersect.cpp"
#endif  // GRID_INTERSECT_H
