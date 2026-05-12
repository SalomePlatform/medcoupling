// Copyright (C) 2007-2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#pragma once

#include "BBTree.txx"
#include "InterpKernelBBoxDistance.txx"

#include <limits>
#include <array>

// fmt: off
template <int dim>
class BBTreeWithBBox
{
   protected:
    std::array<double, 2 * dim> _bbox;

   public:
    const std::array<double, 2 * dim> &getBBox() const { return _bbox; }
};

template <int dim, class ConnType>
class BBTreeClosestCompact : public BBTreeBaseBaseBase<dim, ConnType, BBTreeClosestCompact<dim, ConnType> >,
                             public BBTreeWithBBox<dim>
{
   public:
    BBTreeClosestCompact(int level) : BBTreeBaseBaseBase<dim, ConnType, BBTreeClosestCompact<dim, ConnType> >(level) {}
    // required for walk concept
    bool empty() const { return this->_bbox[0] == std::numeric_limits<double>::max(); }

   public:
    static BBTreeClosestCompact<dim, ConnType> Deserialize(
        const std::vector<bool> &structure, const std::vector<std::array<double, 2 * dim> > &bboxData
    )
    {
        size_t structurePt(0);
        return DeserializeCompactInternal(0, structurePt, structure, bboxData);
    }

   private:
    void fillDeserializeInternal(
        size_t &structurePt,
        const std::vector<bool> &structure,
        const std::vector<std::array<double, 2 * dim> > &bboxData
    )
    {
        this->_bbox = bboxData[structurePt];
        if (!structure[structurePt++])
        {
            return;
        }
        this->_left.reset(new BBTreeClosestCompact<dim, ConnType>(this->_level + 1));
        this->_left->fillDeserializeInternal(structurePt, structure, bboxData);
        this->_right.reset(new BBTreeClosestCompact<dim, ConnType>(this->_level + 1));
        this->_right->fillDeserializeInternal(structurePt, structure, bboxData);
    }

    static BBTreeClosestCompact<dim, ConnType> DeserializeCompactInternal(
        int level,
        size_t &structurePt,
        const std::vector<bool> &structure,
        const std::vector<std::array<double, 2 * dim> > &bboxData
    )
    {
        BBTreeClosestCompact<dim, ConnType> ret(level);
        ret.fillDeserializeInternal(structurePt, structure, bboxData);
        return ret;
    }
};

template <int dim, class ConnType>
void
BBTreeBoundaryCompute(const double *bboxPtr, const std::vector<ConnType> &elems, double boundary[2 * dim])
{
    for (int i = 0; i < dim; ++i)
    {
        boundary[2 * i] = std::numeric_limits<double>::max();
        boundary[2 * i + 1] = -std::numeric_limits<double>::max();
    }
    for (auto elem : elems)
    {
        for (int i = 0; i < dim; ++i)
        {
            boundary[2 * i] = std::min(boundary[2 * i], bboxPtr[elem * 2 * dim + i * 2]);
            boundary[2 * i + 1] = std::max(boundary[2 * i + 1], bboxPtr[elem * 2 * dim + i * 2 + 1]);
        }
    }
}

/*!
 * [EDF34966] Class in charge to find closest elements recursively in \a this relatively to a point.
 * Contrary to bbox intersections, here axis are coupled.
 * The idea is :
 *
 * 0 - split in BBTree without overlaping
 * 1 - compute coarse grain bboxes by aggregation
 * 2 - For a point compute min max for each instance of BBTreeClosest (_min is updated) and deduce the min of max
 * distance (called minOfMaxes). 3 - Keep Only BBTreeClosest instances having _min < minOfMaxes ( first filtering ) 4 -
 * For selected BBTreeClosest instances compute for each element min/max. Store it in _candidates attibute and refine
 * _min attibutes and return a refined minOfMaxes 5 - With refined minOfMaxes at step 4 return elements having min
 * distance <= refined minOfMaxes
 *
 * Steps 2, 3, 4 and 5 are repeated for each point.
 */
template <int dim, class ConnType>
class BBTreeClosest : public BBTreeBaseBase<dim, ConnType, 15, 12, BBTreeClosest<dim, ConnType> >,
                      public BBTreeWithBBox<dim>
{
   private:
    mutable double _min = -std::numeric_limits<double>::max();
    mutable std::vector<std::pair<ConnType, double> > _candidates;

   private:
    void fillDeserializeInternal(
        size_t &structurePt,
        const std::vector<bool> &structure,
        const std::vector<std::array<double, 2 * dim> > &bboxData,
        size_t &elemsDataPt,
        const std::vector<ConnType> &elemsData
    )
    {
        this->_bbox = bboxData[structurePt];
        if (!structure[structurePt++])
        {
            ConnType sz(elemsData[elemsDataPt++]);
            this->_elems.insert(
                this->_elems.end(), elemsData.cbegin() + structurePt, elemsData.cbegin() + structurePt + sz
            );
            elemsDataPt += sz;
            return;
        }
        this->_left.reset(new BBTreeClosest<dim, ConnType>(this->_level + 1));
        this->_left->fillDeserializeInternal(structurePt, structure, bboxData, elemsDataPt, elemsData);
        this->_right.reset(new BBTreeClosest<dim, ConnType>(this->_level + 1));
        this->_right->fillDeserializeInternal(structurePt, structure, bboxData, elemsDataPt, elemsData);
    }

    void fillDeserializeCompactInternal(
        size_t &structurePt,
        const std::vector<bool> &structure,
        const std::vector<std::array<double, 2 * dim> > &bboxData
    )
    {
        this->_bbox = bboxData[structurePt];
        if (!structure[structurePt++])
        {
            return;
        }
        this->_left.reset(new BBTreeClosest<dim, ConnType>(this->_level + 1));
        this->_left->fillDeserializeCompactInternal(structurePt, structure, bboxData);
        this->_right.reset(new BBTreeClosest<dim, ConnType>(this->_level + 1));
        this->_right->fillDeserializeCompactInternal(structurePt, structure, bboxData);
    }

    static BBTreeClosest<dim, ConnType> DeserializeInternal(
        int level,
        size_t &structurePt,
        const std::vector<bool> &structure,
        const std::vector<std::array<double, 2 * dim> > &bboxData,
        size_t &elemsDataPt,
        const std::vector<ConnType> &elemsData
    )
    {
        BBTreeClosest<dim, ConnType> ret(level);
        ret.fillDeserializeInternal(structurePt, structure, bboxData, elemsDataPt, elemsData);
        return ret;
    }

    static std::unique_ptr<double[]> FillClosestDataNotTerminal(
        int level, const double *bbs, const ConnType *elems, ConnType nbelems
    )
    {
        std::unique_ptr<double[]> nodes(new double[nbelems]);
        for (ConnType i = 0; i < nbelems; i++)
        {
            ConnType elem;
            if (elems)
                elem = elems[i];
            else
                elem = i;
            nodes[i] = (bbs[elem * dim * 2 + (level % dim) * 2] + bbs[elem * dim * 2 + (level % dim) * 2 + 1]) / 2.0;
        }
        return nodes;
    }

    static void SplitClosestForNotTerminal(
        int level,
        const double *bbs,
        const ConnType *elems,
        ConnType nbelems,
        std::vector<ConnType> &left,
        std::vector<ConnType> &right
    )
    {
        auto ElemGetter = [elems](ConnType i) -> ConnType
        {
            if (elems)
                return elems[i];
            else
                return i;
        };

        std::unique_ptr<ConnType[]> idx(new ConnType[nbelems]);
        for (ConnType i = 0; i < nbelems; ++i) idx[i] = i;
        std::unique_ptr<double[]> nodes(FillClosestDataNotTerminal(level, bbs, elems, nbelems));
        ConnType median(nbelems / 2);
        std::nth_element(
            idx.get(),
            idx.get() + median,
            idx.get() + nbelems,
            [&](ConnType a, ConnType b) { return nodes[a] < nodes[b]; }
        );
        //
        left.reserve(median);
        right.reserve(nbelems - median);
        for (ConnType i = 0; i < median; ++i) left.push_back(ElemGetter(idx[i]));

        for (ConnType i = median; i < nbelems; ++i) right.push_back(ElemGetter(idx[i]));
    }

    BBTreeClosest(int level) : BBTreeBaseBase<dim, ConnType, 15, 12, BBTreeClosest<dim, ConnType> >(level) {}

   public:
    BBTreeClosest() = default;
    BBTreeClosest(const double *bbs, const ConnType *elems, int level, ConnType nbelems)
        : BBTreeBaseBase<dim, ConnType, 15, 12, BBTreeClosest<dim, ConnType> >(bbs, elems, level, nbelems)
    {
        if (this->constexprTerminal(level, nbelems))
        {
            this->computeBBox();
        }
        else
        {
            std::vector<ConnType> new_elems_left, new_elems_right;
            {
                this->SplitClosestForNotTerminal(level, bbs, elems, nbelems, new_elems_left, new_elems_right);
            }
            auto FromVectToPt = [](const std::vector<ConnType> &v) -> const ConnType *
            { return v.empty() ? nullptr : v.data(); };
            this->_left.reset(new BBTreeClosest<dim, ConnType>(
                bbs, FromVectToPt(new_elems_left), level + 1, (ConnType)new_elems_left.size()
            ));
            this->_right.reset(new BBTreeClosest<dim, ConnType>(
                bbs, FromVectToPt(new_elems_right), level + 1, (ConnType)new_elems_right.size()
            ));

            for (int i = 0; i < dim; ++i)
            {
                this->_bbox[2 * i] = std::min(this->left()->_bbox[2 * i], this->right()->_bbox[2 * i]);
                this->_bbox[2 * i + 1] = std::max(this->left()->_bbox[2 * i + 1], this->right()->_bbox[2 * i + 1]);
            }
        }
    }

    void serialize(
        std::vector<bool> &structure,
        std::vector<std::array<double, 2 * dim> > &bboxData,
        std::vector<ConnType> &elemsData
    ) const
    {
        bboxData.push_back(this->_bbox);
        if (this->terminal())
        {
            structure.push_back(false);
            elemsData.push_back((ConnType)this->_elems.size());
            elemsData.insert(elemsData.end(), this->_elems.cbegin(), this->_elems.cend());
        }
        else
        {
            structure.push_back(true);
            this->_left->serialize(structure, bboxData, elemsData);
            this->_right->serialize(structure, bboxData, elemsData);
        }
    }

    void serializeCompact(std::vector<bool> &structure, std::vector<std::array<double, 2 * dim> > &bboxData) const
    {
        bboxData.push_back(this->_bbox);
        if (this->terminal())
        {
            structure.push_back(false);
        }
        else
        {
            structure.push_back(true);
            this->_left->serializeCompact(structure, bboxData);
            this->_right->serializeCompact(structure, bboxData);
        }
    }

    static BBTreeClosest<dim, ConnType> Deserialize(
        const std::vector<bool> &structure,
        const std::vector<std::array<double, 2 * dim> > &bboxData,
        const std::vector<ConnType> &elemsData
    )
    {
        size_t structurePt(0), elemsDataPt(0);
        return DeserializeInternal(0, structurePt, structure, bboxData, elemsDataPt, elemsData);
    }

    static BBTreeClosest<dim, ConnType> DeserializeCompact(
        const std::vector<bool> &structure, const std::vector<std::array<double, 2 * dim> > &bboxData
    )
    {
        BBTreeClosest<dim, ConnType> ret(0);
        size_t structurePt(0);
        ret.fillDeserializeCompactInternal(structurePt, structure, bboxData);
        return ret;
    }

    /*!
     * Compute Min of Maxes between \a inputBBox and bboxes contained in \a this.
     */
    void bboxMinOfMaxes(const std::array<double, 2 * dim> &inputBBox, double &res) const
    {
        if (this->terminal())
        {
            double zeMin, zeMax;
            INTERP_KERNEL::HighLevelBBDistInternal<dim>(inputBBox, this->_bbox, zeMin, zeMax);
            res = std::min(res, zeMax);
        }
        else
        {
            this->_left->bboxMinOfMaxes(inputBBox, res);
            this->_right->bboxMinOfMaxes(inputBBox, res);
        }
    }

    void bboxSelect(
        const std::array<double, 2 * dim> &inputBBox,
        double thres,
        std::set<const BBTreeClosest<dim, ConnType> *> &blockSelected
    ) const
    {
        for (const auto &leaf : *this)
        {
            const decltype(*this) &leaf2(static_cast<const decltype(*this) &>(leaf));
            double zeMin, zeMax;
            INTERP_KERNEL::HighLevelBBDistInternal<dim>(inputBBox, leaf2.getBBox(), zeMin, zeMax);
            if (zeMin < thres)
            {
                blockSelected.insert(&leaf2);
            }
        }
    }

    /*!
     * Closest candidate between a point and bboxes in this
     */
    void getClosestCandidates(const double pt[dim], std::vector<ConnType> &elems) const
    {
        auto BackInserter = [&elems](mcIdType v) { elems.push_back(v); };
        this->getClosestCandidatesGen(pt, BackInserter);
    }

    template <class T>
    void getClosestCandidatesGen(const double pt[dim], T container) const
    {
        double minOfMaxes(this->getMinOfMaxes(pt));
        minOfMaxes = this->refineMax(pt, minOfMaxes);
        this->filter(minOfMaxes, container);
    }

   public:
    /*!
     * Step 2 of algorithm. Return the minOfMaxes but also update _min attribute
     * non const (_min is updated). First turn. no optimization
     */
    double getMinOfMaxes(const double pt[dim]) const
    {
        if (this->terminal())
        {
            double val;
            INTERP_KERNEL::HighLevelDistInternal<dim>(pt, this->_bbox.data(), _min, val);
            return val;
        }
        else
        {
            double lret(this->left()->getMinOfMaxes(pt));
            double rret(this->right()->getMinOfMaxes(pt));
            _min = std::min(this->left()->_min, this->right()->_min);
            return std::min(lret, rret);
        }
    }
    /*!
     * Step 3 and 4 of algorithm. Refine _min / _max by considering not only the whole bbox of all
     * elements in terminal but each bbox of each element in terminal.
     */
    double refineMax(const double pt[dim], double dist) const
    {
        if (_min > dist)
            return dist;
        if (this->terminal())
        {
            double ret(std::numeric_limits<double>::max());
            const double *bbptr(this->getBBoxData());
            _min = std::numeric_limits<double>::max();
            for (auto elem : this->getElements())
            {
                double zeMin, zeMax;
                INTERP_KERNEL::HighLevelDistInternal<dim>(pt, bbptr + elem * 2 * dim, zeMin, zeMax);
                _min = std::min(_min, zeMin);
                if (zeMin <= dist)
                {
                    _candidates.emplace_back(std::pair<ConnType, double>(elem, zeMin));
                    ret = std::min(ret, zeMax);
                }
            }
            return ret;
        }
        else
        {
            double lval(this->left()->refineMax(pt, dist));
            double rval(this->right()->refineMax(pt, std::min(dist, lval)));  // optimization
            _min = std::min(this->left()->_min, this->right()->_min);
            return std::min(lval, rval);
        }
    }
    /*!
     * Step 5 of algorithm.
     */
    template <class T>
    void filter(double dist, T container) const
    {
        if (_min > dist)
        {
            _candidates.clear();
            return;
        }
        if (this->terminal())
        {
            for (auto elem : this->_candidates)
            {
                if (elem.second <= dist)
                {
                    container(elem.first);
                }
            }
            _candidates.clear();
        }
        else
        {
            this->left()->filter(dist, container);
            this->right()->filter(dist, container);
        }
    }
    void computeBBox()
    {
        const std::vector<ConnType> &elems(this->getElements());
        const double *bboxPtr(this->getBBoxData());
        BBTreeBoundaryCompute<dim, ConnType>(bboxPtr, elems, this->_bbox.data());
    }
};
// fmt: on
